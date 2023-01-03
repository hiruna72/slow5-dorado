#include <metal_stdlib>
using namespace metal;

constant int TILE_SIZE = 8;

// Values set via the FunctionConstantValues object passed in at MTL::Function
// creation time.
constant int kLstmLayerSize [[function_constant(0)]];
constant bool kLstmReversedInTime [[function_constant(1)]];
constant int kLinearContractDim [[function_constant(2)]];
constant int kLinearInnerDim [[function_constant(3)]];
constant bool kConvOutputClamp [[function_constant(4)]];
constant float kLinearOutputScale [[function_constant(5)]];
constant bool kLinearOutputClamp [[function_constant(6)]];
constant bool kLinearOutputTanh [[function_constant(7)]];
constant bool kLinearOutputAsByte [[function_constant(8)]];

namespace {

inline float sigmoid(float x) {
    return 1.f / (1.f + metal::exp(-x));
}

inline float tanh_fast(float x) {
    return 2.f * sigmoid(2.f * x) - 1.f;
}

inline float conv_activation(float x) {
    // SiLU / swish activation.
    const float y = x * sigmoid(x);
    if (kConvOutputClamp) {
        // Note: the lower bound is inoperative, since SiLU(x) has a min. of ~-0.28.
        // Needs to updated anyway, once clamp range is nailed down.
        return clamp(y, -0.5f, 3.5f);
    }
    return y;
}

}

// Precision of input activations and weights (before conversion).
typedef float ftype_in;

// Precision of layer processing.
#if 0
typedef float ftype;
typedef metal::simdgroup_float8x8 simdgroup_ftype8x8;
#else
typedef half ftype;
typedef metal::simdgroup_half8x8 simdgroup_ftype8x8;
#endif

// Precision of back guides and posterior probabilities.
// (Scores are int8_t.)
typedef float ftype_out;

#define MAX_LAYER_SIZE 512
#define KERNEL_INDEX_INPUTS \
    uint tid [[thread_index_in_threadgroup]], \
    uint gid [[threadgroup_position_in_grid]], \
    uint sid [[simdgroup_index_in_threadgroup]], \
    uint simdgroups [[simdgroups_per_threadgroup]], \
    uint threadgroups [[threadgroups_per_grid]], \
    uint threads [[threads_per_threadgroup]]

struct ScanArgs {
    int T;
    int N;
    int C;
};

// Scores must be rescaled from byte range to [-5.0, 5.0] before use in
// forward / backward scans.
float ScaleByteScore(int8_t byte_score) {
    constexpr auto kScoreScale = static_cast<float>(5.0 / 127.0);
    return kScoreScale * static_cast<float>(byte_score);
}

kernel void backward_scan(
    device const ScanArgs* const args,
    device const int8_t* const scores_in,
    device ftype_out* const out,
    KERNEL_INDEX_INPUTS)
{
    constexpr int kNumBases = 4;
    constexpr int kNumTransitions = kNumBases + 1;
    constexpr float kFixedStayScore = 2.0f;

    const int T = args->T;
    const int N = args->N;
    const int num_states = args->C;
    const int ts_states = num_states * kNumBases;
    const int chunk = gid;

    device const int8_t* const chunk_in = scores_in + chunk * ts_states;
    device ftype_out* const chunk_out = out + chunk * (T+1) * num_states;
    device ftype_out* const alpha_init = chunk_out + num_states * T;
    for (int c = tid; c < num_states; c += threads) {
        alpha_init[c] = 0.0f;
    }
    for (int ts = 0; ts < T; ++ts) {
        threadgroup_barrier(mem_flags::mem_device);
        device const auto* const ts_in = chunk_in + N * ts_states * (T - ts - 1);
        device ftype_out* const ts_alpha_in = alpha_init - num_states * ts;
        device ftype_out* const ts_alpha_out = ts_alpha_in - num_states;

        const int state = tid;
        const int stay_state_idx = state;
        const int step_state_idx_a = (state * kNumBases) % num_states;
        const int step_trans_idx_a = step_state_idx_a * kNumBases +
            ((state * kNumBases) / num_states);

        float vals[kNumTransitions];
        float max_val = vals[0] = ts_alpha_in[stay_state_idx] + kFixedStayScore;
        for (int base = 0; base < kNumBases; ++base) {
            vals[base + 1] = ts_alpha_in[step_state_idx_a + base] +
                ScaleByteScore(ts_in[step_trans_idx_a + base * kNumBases]);
            max_val = max(max_val, vals[base + 1]);
        }
        float sum = 0.0f;
        for (int i = 0; i < kNumTransitions; ++i) {
            sum += exp(vals[i] - max_val);
        }
        ts_alpha_out[tid] = max_val + log(sum);
    }
}

kernel void forward_scan(
    device const ScanArgs* const args,
    device const int8_t* const scores_in,
    device ftype_out* const out,
    KERNEL_INDEX_INPUTS)
{
    constexpr int kNumBases = 4;
    constexpr int kNumTransitions = kNumBases + 1;
    constexpr float kFixedStayScore = 2.0f;

    const int T = args->T;
    const int N = args->N;
    const int num_states = args->C;
    const int ts_states = num_states * kNumBases;
    const int kMsb = num_states / kNumBases;
    const int chunk = gid;

    device const int8_t* const chunk_in = scores_in + chunk * ts_states;
    device ftype_out* const chunk_out = out + chunk * (T+1) * num_states;
    device ftype_out* const alpha_init = chunk_out;
    for (int c = tid; c < num_states; c += threads) {
        alpha_init[c] = 0.0f;
    }
    for (int ts = 0; ts < T; ++ts) {
        threadgroup_barrier(mem_flags::mem_device);
        device const auto* const ts_in = chunk_in + N * ts_states * ts;
        device ftype_out* const ts_alpha_in = alpha_init + num_states * ts;
        device ftype_out* const ts_alpha_out = ts_alpha_in + num_states;

        const int state = tid;
        const int stay_state_idx = state;
        const int step_state_idx_a = state / kNumBases;
        const int step_trans_idx_a = state * kNumBases;

        float vals[kNumTransitions];
        float max_val = vals[0] = ts_alpha_in[stay_state_idx] + kFixedStayScore;
        for (int base = 0; base < kNumBases; ++base) {
            vals[base + 1] = ts_alpha_in[step_state_idx_a + base * kMsb] +
                ScaleByteScore(ts_in[step_trans_idx_a + base]);
            max_val = max(max_val, vals[base + 1]);
        }
        float sum = 0.0f;
        for (int i = 0; i < kNumTransitions; ++i) {
            sum += exp(vals[i] - max_val);
        }
        ts_alpha_out[tid] = max_val + log(sum);
    }
}

kernel void add_softmax(
    device const ScanArgs* const args,
    device ftype_out* const fwd_post,
    device const ftype_out* const bwd,
    KERNEL_INDEX_INPUTS)
{
    int T = args->T + 1;
    int C = args->C;
    int chunk = gid;
    int simd_lane = tid & 31;

    for (int ts = sid; ts < T; ts += simdgroups) {
        int ts_idx = (chunk * T + ts) * C;
        float max_val = -1e38;
        for (int i = simd_lane; i < C; i += 32) {
            float val = fwd_post[ts_idx + i] + bwd[ts_idx + i];
            max_val = max(max_val, val);
            fwd_post[ts_idx + i] = val;
        }
        max_val = simd_max(max_val);
        float sum = 0;
        for (int i = simd_lane; i < C; i += 32) {
            float val = exp(fwd_post[ts_idx + i] - max_val);
            sum += val;
            fwd_post[ts_idx + i] = val;
        }
        sum = simd_sum(sum);
        float rcp_sum = 1.f / sum;
        for (int i = simd_lane; i < C; i += 32) {
            fwd_post[ts_idx + i] *= rcp_sum;
        }
    }
}

struct ConvArgs {
    int in_size;
    int win_size;
    int out_size;
    int stride;
    int pad;
    int chunk_size_in; // NOTE: multiple of stride!
    int num_chunks; // Actually batch size
};

/*
// Generic convolution implementation that assumes that `weights` are provided as contiguous tensor of shape
// [W, Cin, Cout] (or [window_size, in_size, out_size]), `in` is a contiguous tensor of shape [N, Tin, Cin] (or
// [num_chunks, chunk_size_in, in_size]), and `out` is a contiguous tensor of shape [N, Tout, Cout] (or
// [num_chunks, chunk_size_in / stride, out_size]).
//
// Should work in place of specialised versions below, provided zero-padding requirements are taken into account, and
// the output of the last convolution layer is passed to `reorder_input` before it is processed by the `lstm` kernel.
kernel void conv(
    device const ConvArgs* const args,
    device const ftype* const in,
    device const ftype* const weights,
    device ftype* const out,
    KERNEL_INDEX_INPUTS)
{
    const int in_size = args->in_size;
    const int win_size = args->win_size;
    const int dp_size = in_size * win_size;
    const int out_size = args->out_size;
    const int stride = args->stride;
    const int pad = args->pad;
    const int chunk_size_in = args->chunk_size_in;
    const int chunk_size_out = chunk_size_in / stride;
    const int num_chunks = args->num_chunks;

    for (int chunk = gid * threads + tid; chunk < num_chunks; chunk += threadgroups * threads) {
        for (int ts = 0; ts < chunk_size_out; ++ts) {
            int in_pos_start = (ts * stride - pad) * in_size;
            for (int output_idx = 0; output_idx < out_size; ++output_idx) {
                ftype sum = weights[dp_size * out_size + output_idx]; // bias
                for (int dp_pos = 0; dp_pos < dp_size; ++dp_pos) {
                    int in_pos = in_pos_start + dp_pos;
                    if (in_pos >= 0 && in_pos < chunk_size_in * in_size) {
                        sum += in[chunk * chunk_size_in * in_size + in_pos] * weights[dp_pos * out_size + output_idx];
                    }
                }
                out[chunk * chunk_size_out * out_size + ts * out_size + output_idx] = conv_activation(sum);
            }
        }
    }
}
*/

// Rearranges/replicates weights for v3-type conv1, where output feature size is 4.
// `weights_in` is expected to be a contiguous tensor of shape [`win_size` + 1, `out_size`],
// i.e. [6, 4], where `weights_in[win_size]` is the bias vector.
// `weights_out` is a contiguous tensor of shape [19, 8]
//
//   auto in = torch::empty({6, 4}, torch::kF32); // ... populated with weights/bias
//   auto out = torch::zeros({19, 8}, torch::kF16);
//   out.index({Slice(6, 11), Slice(0, 4)}) = in.slice(0, 0, 5);
//   out.index({Slice(7, 12), Slice(4, 8)}) = in.slice(0, 0, 5);
//   out.index({18, Slice(0, 4)}) = in[5];
//   out.index({18, Slice(4, 8)}) = in[5];
//
kernel void conv1_out4_simd_reorder_weights(
    device const ConvArgs* const args,
    device const ftype_in* const weights_in,
    device ftype* const weights_out)
{
    const int win_size = 5;
    const int out_size = 4;
    for (int vec4_col = 0; vec4_col < 2; ++vec4_col) {
        for (int in_col = 0; in_col < out_size; ++in_col) {
            for (int row = 0; row < 19; ++row) {
                int in_row = (row >= 11 + vec4_col) ? ((row == 18) ? win_size : -1) : row - 6 - vec4_col;
                weights_out[row * TILE_SIZE + vec4_col * 4 + in_col] =
                    (in_row >= 0) ? weights_in[in_row * out_size + in_col] : ftype(0);
            }
        }
    }
}

// Rearranges/replicates weights for a conv1 where output feature size is a multiple of TILE_SIZE.
// `weights_in` is expected to be a contiguous tensor of shape [`win_size` + 1, `out_size`],
// where `weights_in[win_size]` is the bias vector.
// `weights_out` is a contiguous tensor of shape [`win_size` + 2 * TILE_SIZE - 1, `out_size`]
//
//   auto in = torch::empty({win_size + 1, out_size}, torch::kF32); // ... populated with weights/bias
//   auto out = torch::zeros({win_size + 2 * TILE_SIZE - 1, out_size}, torch::kF16);
//   out.slice(0, TILE_SIZE - 1, TILE_SIZE + win_size - 1) = in.slice(0, 0, win_size);
//   out[win_size + 2 * TILE_SIZE - 2] = in[win_size];
//
kernel void conv1_out16_simd_reorder_weights(
    device const ConvArgs* const args,
    device const ftype_in* const weights_in,
    device ftype* const weights_out)
{
    const int win_size = args->win_size;
    const int out_size = args->out_size;
    const int rows = win_size + 2 * TILE_SIZE - 1;
    for (int col = 0; col < out_size; ++col) {
        for (int row = 0; row < rows; ++row) {
            int in_row = (row >= TILE_SIZE + win_size - 1) ? ((row == rows - 1) ? win_size : -1) : row - TILE_SIZE + 1;
            weights_out[row * out_size + col] = (in_row >= 0) ? weights_in[in_row * out_size + col] : ftype(0);
        }
    }
}

// Rearranges/replicates weights for v3-type conv2, where input feature size is 4.
//
//   auto in = torch::empty({21, 16}, torch::kF32); // ... populated with weights/bias
//   auto out = torch::zeros({29, 16}, torch::kF16);
//   out.slice(0, 4, 24) = in.slice(0, 0, 20);
//   out[28] = in[20];
//
kernel void conv2_in4_simd_reorder_weights(
    device const ConvArgs* const args,
    device const ftype_in* const weights_in,
    device ftype* const weights_out)
{
    for (int col = 0; col < 16; ++col) {
        for (int row = 0; row < 29; ++row) {
            int in_row = (row >= 24) ? ((row == 28) ? 20 : -1) : row - 4;
            weights_out[row * 16 + col] = (in_row >= 0) ? weights_in[in_row * 16 + col] : ftype(0);
        }
    }
}

// Rearranges/replicates weights for v4-type conv2, where input feature size is a multiple of TILE_SIZE.
// Just type conversion
kernel void conv2_in16_simd_reorder_weights(
    device const ConvArgs* const args,
    device const ftype_in* const weights_in,
    device ftype* const weights_out)
{
    const int cols = args->out_size;
    const int rows = args->in_size * args->win_size + 1;
    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < cols; ++col) {
            weights_out[row * cols + col] = weights_in[row * cols + col];
        }
    }
}

// Just type conversion
kernel void conv3_simd_reorder_weights(
    device const ConvArgs* const args,
    device const ftype_in* const weights_in,
    device ftype* const weights_out)
{
    const int cols = args->out_size;
    const int rows = args->in_size * args->win_size + 1;
    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < cols; ++col) {
            weights_out[row * cols + col] = weights_in[row * cols + col];
        }
    }
}

// Just type conversion
kernel void float_to_half(
    device const int* const num_elems,
    device const float* const in,
    device half* const out,
    KERNEL_INDEX_INPUTS)
{
    for (int elem = gid * threads + tid; elem < *num_elems; elem += threadgroups * threads) {
        out[elem] = in[elem];
    }
}


template <int SIMD_TILES_M, int SIMD_TILES_N, int transpose_A> struct MatMul {
    static_assert(SIMD_TILES_M <= 6, "SIMD_TILES_M must be <= 6");
    simdgroup_ftype8x8 A, B[SIMD_TILES_N], accum[SIMD_TILES_M][SIMD_TILES_N];

    void mm_bias(
            int k_tiles_begin, int k_tiles_end,
            device ftype const *a_ptr, int a_stride, int a_col, int a_row,
            device ftype const *b_ptr, int b_stride, int b_col, int b_row,
            device ftype const *bias)
    {
        a_ptr += a_row * a_stride + a_col;
        b_ptr += b_row * b_stride + b_col;
        for (int i = 0; i < SIMD_TILES_N; ++i) {
            for (int j = 0; j < SIMD_TILES_M; ++j) {
                simdgroup_load(accum[j][i], bias + b_col + i * TILE_SIZE, 0);
            }
        }
        for (int k_tile = k_tiles_begin; k_tile < k_tiles_end; ++k_tile) {
            for (int i = 0; i < SIMD_TILES_N; ++i) {
                simdgroup_load(B[i], b_ptr, b_stride,
                               ulong2(i * TILE_SIZE, k_tile * TILE_SIZE));
            }
#define SMAC_ROW(m_tile)\
    if (m_tile < SIMD_TILES_M) {\
        simdgroup_load(A, a_ptr, a_stride, transpose_A ? ulong2(m_tile * TILE_SIZE, k_tile * TILE_SIZE) : ulong2(k_tile * TILE_SIZE, m_tile * TILE_SIZE)); \
        for (int n_tile = 0; n_tile < SIMD_TILES_N; ++n_tile) {\
            simdgroup_multiply_accumulate(accum[m_tile][n_tile], A, B[n_tile], accum[m_tile][n_tile]);\
        }\
    }
            SMAC_ROW(0);
            SMAC_ROW(1);
            SMAC_ROW(2);
            SMAC_ROW(3);
            SMAC_ROW(4);
            SMAC_ROW(5);
#undef SMAC_ROW
        }
    }
};


// Apply conv_activation to a simdgroup matrix tile (using threadgroup memory) then store the tile in device memory.
void conv_activation_and_store(simdgroup_ftype8x8 A, threadgroup ftype *simd_local_buf, int tid,
    device ftype* const out_buf, int out_stride, int tile_col, int tile_row)
{
    simdgroup_store(A, simd_local_buf, TILE_SIZE);
    for (int elem = tid & 31; elem < TILE_SIZE * TILE_SIZE; elem += 32) {
        ftype val = simd_local_buf[elem];
        simd_local_buf[elem] = conv_activation(val);
    }
    simdgroup_load(A, simd_local_buf, TILE_SIZE);
    simdgroup_store(A, out_buf, out_stride, ulong2(tile_col * TILE_SIZE, tile_row * TILE_SIZE));

}


// Specialised conv1 implementation for v3-type models, where output feature size is 4.
// Given a contiguous input tensor of shape [batch_size, chunk_size, 1] this will fill
// a contiguous output tensor of shape [batch_size, chunk_size + 8, 4], where the actual output
// is located at output.slice(1, 2, chunk_size + 2) and values before and after that slice are
// set to zero (this padding is used by the second layer in order to avoid special handling
// of the edges).
#define SIMD_GROUPS 16
[[max_total_threads_per_threadgroup(SIMD_GROUPS * 32)]]
kernel void conv1_out4_simd(
    device const ConvArgs* const args,
    device const ftype* const in_buf,
    device const ftype* const weights_buf,
    device ftype* const out_buf,
    KERNEL_INDEX_INPUTS)
{
//    const int in_size = 1;
    const int out_size = 4;
    const int chunk_size = args->chunk_size_in; // must be multiple of 8
    const int chunk_tiles = args->num_chunks / TILE_SIZE; // num_chunks must be multiple of TILE_SIZE
    const int out_stride = (chunk_size + 8) * out_size;
    threadgroup ftype simd_out_buf[SIMD_GROUPS][TILE_SIZE * TILE_SIZE];
    simdgroup_ftype8x8 W[6], B, I[2], A[4];

    for (int i = 0; i < 6; ++i) {
        simdgroup_load(W[i], weights_buf, TILE_SIZE, ulong2(0, 10 - 2 * i));
    }
    simdgroup_load(B, weights_buf + 18 * TILE_SIZE, 0);
    int num_iters = (chunk_size / TILE_SIZE) - 1;

    // Deal with the chunk edges first
    for (int tile_row = gid; tile_row < chunk_tiles; tile_row += threadgroups) {
        if (sid < 2) {
            int is_last = sid;
            simdgroup_load(I[0], in_buf, chunk_size, ulong2(is_last * (chunk_size - 8), tile_row * TILE_SIZE));
            if (!is_last) {
                // Start of time span / output feature row.
                // Padded with 1 tile = 8 entries.
                A[0] = simdgroup_ftype8x8(0);
                simdgroup_multiply_accumulate(A[1], I[0], W[1], B);
                simdgroup_multiply_accumulate(A[2], I[0], W[2], B);
                simdgroup_multiply_accumulate(A[3], I[0], W[3], B);
            } else {
                // End of time span / output feature row.
                // Padded with 3 tiles = 8 entries.
                simdgroup_multiply_accumulate(A[0], I[0], W[4], B);
                A[1] = simdgroup_ftype8x8(0);
                A[2] = simdgroup_ftype8x8(0);
                A[3] = simdgroup_ftype8x8(0);
            }
            for (int i = 0; i < 4; ++i) {
                conv_activation_and_store(A[i], simd_out_buf[sid], tid, out_buf, out_stride, is_last * (num_iters + 1) * 4 + i, tile_row);
            }
        }
    }

    for (int tile_row = gid; tile_row < chunk_tiles; tile_row += threadgroups) {
        for (int iter = sid; iter < num_iters; iter += simdgroups) {
            simdgroup_load(I[0], in_buf, chunk_size, ulong2(iter * TILE_SIZE, tile_row * TILE_SIZE));
            simdgroup_load(I[1], in_buf, chunk_size, ulong2((iter + 1) * TILE_SIZE, tile_row * TILE_SIZE));
            simdgroup_multiply_accumulate(A[0], I[0], W[4], B);
            simdgroup_multiply_accumulate(A[1], I[0], W[5], B);
            simdgroup_multiply_accumulate(A[2], I[1], W[2], B);
            simdgroup_multiply_accumulate(A[3], I[1], W[3], B);
            simdgroup_multiply_accumulate(A[0], I[1], W[0], A[0]);
            simdgroup_multiply_accumulate(A[1], I[1], W[1], A[1]);
            for (int i = 0; i < 4; ++i) {
                conv_activation_and_store(A[i], simd_out_buf[sid], tid, out_buf, out_stride, (iter + 1) * 4 + i, tile_row);
            }
        }
    }
}

#undef SIMD_GROUPS

// Conv1 implementation for v4-type models, where output feature size is 16.
#define SIMD_GROUPS 16
[[max_total_threads_per_threadgroup(SIMD_GROUPS * 32)]]
kernel void conv1_out16_simd(
    device const ConvArgs* const args,
    device const ftype* const in_buf,
    device const ftype* const weights_buf,
    device ftype* const out_buf,
    KERNEL_INDEX_INPUTS)
{
//    const int in_size = 1;
    const int out_size = 16;
    const int chunk_size = args->chunk_size_in; // must be multiple of 8
    const int chunk_tiles = args->num_chunks / TILE_SIZE; // num_chunks must be multiple of TILE_SIZE
    const int out_stride = chunk_size * out_size;
    threadgroup ftype simd_out_buf[SIMD_GROUPS][TILE_SIZE * TILE_SIZE];
    simdgroup_ftype8x8 W[2][12], B[2], I[2], A[8];

    for (int tile_col = 0; tile_col < 2; ++tile_col) {
        for (int tile_row = 0; tile_row < 12; ++tile_row) {
            simdgroup_load(W[tile_col][tile_row], weights_buf, 2 * TILE_SIZE, ulong2(tile_col * TILE_SIZE, tile_row));
        }
        simdgroup_load(B[tile_col], weights_buf + 19 * out_size + tile_col * TILE_SIZE, 0);
    }
    int num_iters = (chunk_size / TILE_SIZE) - 1;

    // Deal with the chunk edges first
    for (int tile_row = gid; tile_row < chunk_tiles; tile_row += threadgroups) {
        if (sid < 2) {
            simdgroup_load(I[0], in_buf, chunk_size, ulong2(0, tile_row * TILE_SIZE));
            simdgroup_load(I[1], in_buf, chunk_size, ulong2(chunk_size - TILE_SIZE, tile_row * TILE_SIZE));
            simdgroup_multiply_accumulate(A[0], I[0], W[sid][9], B[sid]);
            simdgroup_multiply_accumulate(A[1], I[0], W[sid][8], B[sid]);
            simdgroup_multiply_accumulate(A[2], I[0], W[sid][7], B[sid]);
            simdgroup_multiply_accumulate(A[3], I[0], W[sid][6], B[sid]);
            simdgroup_multiply_accumulate(A[4], I[0], W[sid][5], B[sid]);
            simdgroup_multiply_accumulate(A[5], I[0], W[sid][4], B[sid]);
            simdgroup_multiply_accumulate(A[6], I[1], W[sid][3], B[sid]);
            simdgroup_multiply_accumulate(A[7], I[1], W[sid][2], B[sid]);
            for (int i = 0; i < 8; ++i) {
                int tile_col = (int(i > 5) * (chunk_size - 8) + i) * (out_size / TILE_SIZE) + sid;
                conv_activation_and_store(A[i], simd_out_buf[sid], tid, out_buf, out_stride, tile_col, tile_row);
            }
        }
    }

    for (int tile_row = gid; tile_row < chunk_tiles; tile_row += threadgroups) {
        for (int iter = sid; iter < num_iters; iter += simdgroups) {
            simdgroup_load(I[0], in_buf, chunk_size, ulong2(iter * TILE_SIZE, tile_row * TILE_SIZE));
            simdgroup_load(I[1], in_buf, chunk_size, ulong2((iter + 1) * TILE_SIZE, tile_row * TILE_SIZE));
            for (int i = 0; i < 2; ++i) {
                simdgroup_multiply_accumulate(A[0], I[1], W[i][11], B[i]);
                simdgroup_multiply_accumulate(A[1], I[1], W[i][10], B[i]);
                simdgroup_multiply_accumulate(A[2], I[1], W[i][9], B[i]);
                simdgroup_multiply_accumulate(A[3], I[1], W[i][8], B[i]);
                simdgroup_multiply_accumulate(A[4], I[1], W[i][7], B[i]);
                simdgroup_multiply_accumulate(A[5], I[1], W[i][6], B[i]);
                simdgroup_multiply_accumulate(A[6], I[1], W[i][5], B[i]);
                simdgroup_multiply_accumulate(A[7], I[1], W[i][4], B[i]);
                simdgroup_multiply_accumulate(A[0], I[0], W[i][3], A[0]);
                simdgroup_multiply_accumulate(A[1], I[0], W[i][2], A[1]);
                simdgroup_multiply_accumulate(A[2], I[0], W[i][1], A[2]);
                simdgroup_multiply_accumulate(A[3], I[0], W[i][0], A[3]);
                for (int j = 0; j < 8; ++j) {
                    int tile_col = (iter * 8 + 6 + j) * (out_size / TILE_SIZE) + i;
                    conv_activation_and_store(A[j], simd_out_buf[sid], tid, out_buf, out_stride, tile_col, tile_row);
                }
            }
        }
    }
}
#undef SIMD_GROUPS

// Specialised conv2 implementation for v3-type models, where input feature size is 4.
#define SIMD_GROUPS 16
[[max_total_threads_per_threadgroup(SIMD_GROUPS * 32)]]
kernel void conv2_in4_simd
(
    device const ConvArgs* const args,
    device const ftype* const in_buf,
    device const ftype* const weights_buf,
    device ftype* const out_buf,
    KERNEL_INDEX_INPUTS
) {
    const int in_size = 4;
    const int out_size = 16;
    const int chunk_size = args->chunk_size_in; // must be multiple of 2
    const int chunk_tiles = args->num_chunks / TILE_SIZE; // num_chunks must be multiple of TILE_SIZE
    const int in_stride = (chunk_size + 8) * in_size;
    const int out_stride = chunk_size * out_size;
    threadgroup ftype simd_out_buf[SIMD_GROUPS][TILE_SIZE * TILE_SIZE];
    simdgroup_ftype8x8 W[3][4], B[2], I[3], A[4];
    device const ftype* b = weights_buf + 28 * 16;

    simdgroup_load(B[0], b, 0);
    simdgroup_load(B[1], b + TILE_SIZE, 0);
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 3; ++j) {
            simdgroup_load(W[j][i * 2 + 0], weights_buf, 2 * TILE_SIZE, ulong2(0 * TILE_SIZE, j * TILE_SIZE + (1 - i) * 4));
            simdgroup_load(W[j][i * 2 + 1], weights_buf, 2 * TILE_SIZE, ulong2(1 * TILE_SIZE, j * TILE_SIZE + (1 - i) * 4));
        }
    }

    for (int tile_row = gid; tile_row < chunk_tiles; tile_row += threadgroups) {
        for (int iter = sid; iter < chunk_size / 2; iter += simdgroups) {
            simdgroup_load(I[0], in_buf, in_stride, ulong2((iter + 0) * TILE_SIZE, tile_row * TILE_SIZE));
            simdgroup_load(I[1], in_buf, in_stride, ulong2((iter + 1) * TILE_SIZE, tile_row * TILE_SIZE));
            simdgroup_load(I[2], in_buf, in_stride, ulong2((iter + 2) * TILE_SIZE, tile_row * TILE_SIZE));
            for (int i = 0; i < 4; ++i) {
                simdgroup_multiply_accumulate(A[i], I[0], W[0][i], B[i & 1]);
                simdgroup_multiply_accumulate(A[i], I[1], W[1][i], A[i]);
                simdgroup_multiply_accumulate(A[i], I[2], W[2][i], A[i]);
            }
            for (int i = 0; i < 4; ++i) {
                conv_activation_and_store(A[i], simd_out_buf[sid], tid, out_buf, out_stride, iter * 4 + i, tile_row);
            }
        }
    }
}
#undef SIMD_GROUPS

#define SIMD_TILES_M 6
#define SIMD_TILES_N 2
#define SIMD_GROUPS 4

[[max_total_threads_per_threadgroup(SIMD_GROUPS * 32)]]
kernel void conv2_in16_simd
(
    device const ConvArgs* const args,
    device const ftype* const in_buf,
    device const ftype* const weights_buf,
    device ftype* const out_buf,
    KERNEL_INDEX_INPUTS
) {
    const int in_size = args->in_size;
    const int win_size = args->win_size;
    const int dp_size = in_size * win_size;
    const int out_size = 16; // required! //args->out_size;
    const int stride = args->stride;
    const int pad = args->pad;
    const int w_pad_rows = 0;
    const int chunk_size_in = args->chunk_size_in;
    const int chunk_size_out = chunk_size_in / stride;
    const int num_chunks = args->num_chunks;
    const int m_blks = num_chunks / (TILE_SIZE * SIMD_TILES_M);
    const int k_tiles = (dp_size + TILE_SIZE - 1) / TILE_SIZE;
    threadgroup ftype simd_out_buf[SIMD_GROUPS][SIMD_TILES_N * SIMD_TILES_M][TILE_SIZE * TILE_SIZE];
    device const ftype* bias = weights_buf + (dp_size + 2 * w_pad_rows) * out_size;
    const int in_buf_stride = chunk_size_in * in_size;
    const int out_buf_stride = chunk_size_out * out_size;
    MatMul<SIMD_TILES_M, SIMD_TILES_N, 0> mm;

    for (int m_blk = gid; m_blk < m_blks; m_blk += threadgroups) {
        for (int ts = sid; ts < chunk_size_in / stride; ts += simdgroups) {
            int start_pos = (ts * stride - pad) * in_size;
            int end_pos = start_pos + dp_size;
            int clamped_start_pos = max(0, start_pos);
            int clamped_end_pos = min(chunk_size_in * in_size, end_pos);
            int start_pad = clamped_start_pos - start_pos;
            int end_pad = end_pos - clamped_end_pos;
            int start_pad_tiles = start_pad / TILE_SIZE;
            start_pad -= start_pad_tiles * TILE_SIZE;
            int w_row_offset = w_pad_rows - start_pad;
            int end_pad_tiles = end_pad / TILE_SIZE;
            int pad_tiles = start_pad_tiles + end_pad_tiles;
            mm.mm_bias(0, k_tiles - pad_tiles,
                in_buf, in_buf_stride, clamped_start_pos, m_blk * SIMD_TILES_M * TILE_SIZE,
                weights_buf, out_size, 0, w_row_offset, bias);
            for (int i = 0; i < SIMD_TILES_M; ++i) {
                for (int j = 0; j < SIMD_TILES_N; ++j) {
                    simdgroup_store(mm.accum[i][j], simd_out_buf[sid][i * SIMD_TILES_N + j], TILE_SIZE);
                }
            }
            simdgroup_barrier(mem_flags::mem_threadgroup);
            for (uint elem = tid & 31; elem < SIMD_TILES_N * SIMD_TILES_M * TILE_SIZE * TILE_SIZE; elem += 32) {
                // swish activation
                ftype val = simd_out_buf[sid][0][elem];
                simd_out_buf[sid][0][elem] = conv_activation(val);
            }
            simdgroup_barrier(mem_flags::mem_threadgroup);
            for (int i = 0; i < SIMD_TILES_M; ++i) {
                int out_row = (m_blk * SIMD_TILES_M + i) * TILE_SIZE;
                for (int j = 0; j < SIMD_TILES_N; ++j) {
                    int tile_idx = i * SIMD_TILES_N + j;
                    int out_col = ts * out_size + j * TILE_SIZE;
                    simdgroup_load(mm.A, simd_out_buf[sid][tile_idx], TILE_SIZE);
                    simdgroup_store(mm.A, out_buf, out_buf_stride, ulong2(out_col, out_row));
                }
            }
        }
    }
}
#undef SIMD_GROUPS
#undef SIMD_TILES_M
#undef SIMD_TILES_N

#define SIMD_TILES_M 6
#define SIMD_TILES_N 4
#define SIMD_GROUPS 4

[[max_total_threads_per_threadgroup(SIMD_GROUPS * 32)]]
kernel void conv3_simd
(
    device const ConvArgs* const args,
    device const ftype* const in_buf,
    device const ftype* const weights_buf,
    device ftype* const out_buf,
    KERNEL_INDEX_INPUTS
) {
    const int in_size = args->in_size;
    const int in_size_tiles = args->in_size / TILE_SIZE;
    const int win_size = args->win_size;
    const int dp_size = in_size * win_size;
    const int out_size = args->out_size;
    const int stride = args->stride;
    const int pad = args->pad;
    const int chunk_size_in = args->chunk_size_in;
    const int chunk_size_out = chunk_size_in / stride;
    const int num_chunks = args->num_chunks;
    const int m_blks = num_chunks / (TILE_SIZE * SIMD_TILES_M);
    const int n_blks = out_size / (TILE_SIZE * SIMD_TILES_N);
    const int k_blks = dp_size / TILE_SIZE;
    threadgroup ftype simd_out_buf[SIMD_GROUPS][SIMD_TILES_N * SIMD_TILES_M][TILE_SIZE * TILE_SIZE];
    device const ftype* bias = weights_buf + dp_size * out_size;
    const int in_buf_stride = chunk_size_in * in_size;
    MatMul<SIMD_TILES_M, SIMD_TILES_N, 0> mm;

    for (int m_blk = gid; m_blk < m_blks; m_blk += threadgroups) {
        for (int chunk = tid; chunk < SIMD_TILES_M * TILE_SIZE; chunk += threads) {
            for (int i = 0; i < out_size; ++i) {
                int idx = i * num_chunks + m_blk * SIMD_TILES_M * TILE_SIZE + chunk;
                out_buf[idx] = 0;
                out_buf[idx + chunk_size_out * out_size] = 0;
            }
        }
    }

    for (int m_blk = gid; m_blk < m_blks; m_blk += threadgroups) {
        for (int ts = 0; ts < chunk_size_out; ++ts) {
            int start_pos = ts * stride - pad;
            int start_pad_tiles = max(0, -start_pos) * in_size_tiles;
            int end_pad_tiles = max(0, start_pos + win_size - chunk_size_in) * in_size_tiles;
            device ftype* out = out_buf + (ts + 1) * num_chunks * out_size; // One timestep of padding as required by LSTM
            for (int n_blk = sid; n_blk < n_blks; n_blk += simdgroups) {
                mm.mm_bias(start_pad_tiles, k_blks - end_pad_tiles,
                    in_buf, in_buf_stride, start_pos * in_size_tiles * TILE_SIZE, m_blk * SIMD_TILES_M * TILE_SIZE,
                    weights_buf, out_size, n_blk * SIMD_TILES_N * TILE_SIZE, 0, bias);
                for (int i = 0; i < SIMD_TILES_M; ++i) {
                    for (int j = 0; j < SIMD_TILES_N; ++j) {
                        simdgroup_store(mm.accum[i][j], simd_out_buf[sid][i * SIMD_TILES_N + j], TILE_SIZE);
                    }
                }
                simdgroup_barrier(mem_flags::mem_threadgroup);
                for (uint elem = tid & 31; elem < SIMD_TILES_N * SIMD_TILES_M * TILE_SIZE * TILE_SIZE; elem += 32) {
                    // swish activation
                    ftype val = simd_out_buf[sid][0][elem];
                    simd_out_buf[sid][0][elem] = conv_activation(val);
                }
                simdgroup_barrier(mem_flags::mem_threadgroup);
                for (int i = 0; i < SIMD_TILES_M; ++i) {
                    uint out_col = (m_blk * SIMD_TILES_M + i) * TILE_SIZE;
                    for (int j = 0; j < SIMD_TILES_N; ++j) {
                        int tile_idx = i * SIMD_TILES_N + j;
                        uint out_row = (n_blk * SIMD_TILES_N + j) * TILE_SIZE;
                        simdgroup_load(mm.A, simd_out_buf[sid][tile_idx], TILE_SIZE);
                        simdgroup_store(mm.A, out, num_chunks, ulong2(out_col, out_row));
                    }
                }
            }
        }
    }
}
#undef SIMD_GROUPS



// auto W = torch::empty({layer_size, 4, layer_size}, torch::kF32);  // {in_size, gate, out_size}, gate order: GIFO
// auto U = torch::empty({layer_size, 4, layer_size}, torch::kF32);  // as above
// auto b = torch::empty({4, layer_size}, torch::kF32);              // {gate, out_size}
// [ ... fill W, U, and b ... ]
// auto weights_buf = torch::empty{2 * layer_size + 1, layer_size, 4}, torch::kF16);
// weights_buf.slice(0, 0, layer_size) = U.transpose(1, 2);
// weights_buf.slice(0, layer_size, 2 * layer_size) = W.transpose(1, 2);
// weights_buf[2 * layer_size] = b.transpose(0, 1);
//
kernel void reorder_lstm_weights(
    device const ftype_in* const W,
    device const ftype_in* const U,
    device const ftype_in* const b,
    device ftype* const weights_buf,
    KERNEL_INDEX_INPUTS)
{
    const int stride = kLstmLayerSize * 4;
    device const ftype_in* const inputs[3] = { U, W, b };
    const int in_rows[3] = { kLstmLayerSize, kLstmLayerSize, 1 };

    for (int m = 0; m < 3; ++m) {
        for (int r = gid; r < in_rows[m]; r += threadgroups) {
            for (int c = tid; c < kLstmLayerSize; c += threads) {
                for (int gate = 0; gate < 4; ++gate) {
                    weights_buf[(m * kLstmLayerSize + r) * stride + c * 4 + gate] = ftype(inputs[m][r * stride + gate * kLstmLayerSize + c]);
                }
            }
        }
    }
}

struct LstmArgs {
    int batch_tiles;
    int chunk_size;
};

// This takes a contiguous input tensor of shape [T, N, C], and contiguous output tensor
// of shape [T+2, C/8, 8, N/8, 8], and performs these assignments:
//     output.slice(0, 1, T+1) = input.view({T, N/8, 8, C/8, 8}.transpose(1, 3);
//     output[0] = 0;
//     output[T+1] = 0;
kernel void reorder_input(
    device const LstmArgs* const args,
    device const ftype_in* const in,
    device ftype* const out,
    KERNEL_INDEX_INPUTS)
{
    threadgroup ftype bfr[MAX_LAYER_SIZE * TILE_SIZE];
    const int layer_tiles = kLstmLayerSize / TILE_SIZE;
    const int batch_tiles = args->batch_tiles;
    const int chunk_size = args->chunk_size;
    for (int batch_tile = gid; batch_tile < batch_tiles; batch_tile += threadgroups) {
        // note: at timestep=-1 and timestep=chunk_size we do zero-padding in order to avoid having to deal with the edges differently
        for (int timestep = -1; timestep <= chunk_size; ++timestep) {
            threadgroup_barrier(mem_flags::mem_threadgroup);
            for (int chunk = 0; chunk < TILE_SIZE; ++chunk) {
                for (int col = tid; col < kLstmLayerSize; col += threads) {
                    const int idx = (timestep * batch_tiles * TILE_SIZE + (batch_tile * TILE_SIZE + chunk)) * kLstmLayerSize + col;
                    const ftype val = (timestep >= 0 && timestep < chunk_size) ? ftype(in[idx]) : ftype(0);
                    bfr[chunk * MAX_LAYER_SIZE + col] = val;
                }
            }
            threadgroup_barrier(mem_flags::mem_threadgroup);
            for (int tile = sid; tile < layer_tiles; tile += simdgroups) {
                simdgroup_ftype8x8 A;
                simdgroup_load(A, bfr + tile * TILE_SIZE, MAX_LAYER_SIZE);
                simdgroup_store(A, out, batch_tiles * TILE_SIZE, ulong2(batch_tile * TILE_SIZE, (timestep + 1) * kLstmLayerSize + tile * TILE_SIZE));
            }
        }
    }
}

// This takes a contiguous input tensor of shape [T+2, C/8, 8, N/8, 8], and contiguous output tensor
// of shape [T, N, C], and performs this copy:
//     output.view({T, N/8, 8, C/8, 8}) = input.slice(0, 1, T+1).transpose(1, 3);
// Note that it ignores input[0,:,:,:] and input[T+1,:,:,:]
kernel void reorder_output(
    device const LstmArgs* const args,
    device const ftype* const in,
    device ftype_out* const out,
    KERNEL_INDEX_INPUTS)
{
    threadgroup ftype bfr[MAX_LAYER_SIZE * TILE_SIZE];
    const int layer_tiles = kLstmLayerSize / TILE_SIZE;
    const int batch_tiles = args->batch_tiles;
    const int chunk_size = args->chunk_size;
    for (int batch_tile = gid; batch_tile < batch_tiles; batch_tile += threadgroups) {
        for (int timestep = 0; timestep < chunk_size; ++timestep) {
            threadgroup_barrier(mem_flags::mem_threadgroup);
            for (int tile = sid; tile < layer_tiles; tile += simdgroups) {
                simdgroup_ftype8x8 A;
                simdgroup_load(A, in, batch_tiles * TILE_SIZE, ulong2(batch_tile * TILE_SIZE, (timestep + 1) * kLstmLayerSize + tile * TILE_SIZE));
                simdgroup_store(A, bfr + tile * TILE_SIZE, MAX_LAYER_SIZE);
            }
            threadgroup_barrier(mem_flags::mem_threadgroup);
            for (int chunk = 0; chunk < TILE_SIZE; ++chunk) {
                for (int col = tid; col < kLstmLayerSize; col += threads) {
                    const int idx = (timestep * batch_tiles * TILE_SIZE + (batch_tile * TILE_SIZE + chunk)) * kLstmLayerSize + col;
                    out[idx] = ftype_out(bfr[chunk * MAX_LAYER_SIZE + col]);
                }
            }
        }
    }
}

// Note: max_total_threads_per_threadgroup is set via ComputePipelineDescriptor,
// rather than an attribute here, since it depends on the SIMD group count,
// which varies according to LSTM layer size.
kernel void lstm(
        device const LstmArgs* const args,
        device ftype* const in_out,
        device const ftype* const weights_buf,
        device ftype* const state_buf,
        device ftype* const temp_result_buf,
        // The sizes of these buffers are set via MTL::ComputeCommandEncoder.
        // They depend on the SIMD group count.
        threadgroup ftype (* const simd_res_buf)[2 * TILE_SIZE * TILE_SIZE],
        threadgroup ftype (* const simd_out_buf)[TILE_SIZE * TILE_SIZE],
        KERNEL_INDEX_INPUTS)
{
    const int chunk_size = args->chunk_size;
    const int batch_tiles = args->batch_tiles;
    const int m_blks = batch_tiles / SIMD_TILES_M;
    const int n_blks = kLstmLayerSize * 4 / (TILE_SIZE * SIMD_TILES_N);
    const int k_tiles = kLstmLayerSize * 2 / TILE_SIZE;
    const int inout_stride = batch_tiles * TILE_SIZE;
    const int w_stride = kLstmLayerSize * 4;
    MatMul<SIMD_TILES_M, SIMD_TILES_N, 1> mm;
    device const ftype* const bias = weights_buf + 2 * kLstmLayerSize * w_stride;

    const uint t_idx = tid & 31;
    const uint col_bits = t_idx & 3;
    const uint row = t_idx >> 2;
    const uint rb_idx = t_idx * 4;

    for (int m_blk = gid; m_blk < m_blks; m_blk += threadgroups) {
        for (int chunk = tid; chunk < SIMD_TILES_M * TILE_SIZE; chunk += threads) {
            for (int i = 0; i < kLstmLayerSize; ++i) {
                state_buf[i * batch_tiles * TILE_SIZE + m_blk * SIMD_TILES_M * TILE_SIZE + chunk] = 0;
            }
        }
    }

    for (int iter = 0; iter < chunk_size; ++iter) {
        threadgroup_barrier(mem_flags::mem_device | mem_flags::mem_threadgroup);
        const int timestep_out = kLstmReversedInTime ? chunk_size - iter : iter + 1;
        const int timestep_in = kLstmReversedInTime ? timestep_out : timestep_out - 1;
        device const ftype* const in = in_out + timestep_in * inout_stride * kLstmLayerSize;
        device ftype* const out = in_out + timestep_out * inout_stride * kLstmLayerSize;
        for (int m_blk = gid; m_blk < m_blks; m_blk += threadgroups) {
            for (int n_blk = sid; n_blk < n_blks; n_blk += simdgroups) {
                mm.mm_bias(0, k_tiles, in, inout_stride, m_blk * SIMD_TILES_M * TILE_SIZE, 0,
                           weights_buf, w_stride, n_blk * SIMD_TILES_N * TILE_SIZE, 0, bias);
                for (int i = 0; i < SIMD_TILES_M; ++i) {
                    const uint out_chunk_base = (m_blk * SIMD_TILES_M + i) * TILE_SIZE;
                    const uint chunk_idx = out_chunk_base + row;
                    for (int j = 0; j < SIMD_TILES_N; j += 2) {
                        simdgroup_store(mm.accum[i][j + 0], simd_res_buf[sid],
                                        2 * TILE_SIZE);
                        simdgroup_store(mm.accum[i][j + 1], simd_res_buf[sid] + TILE_SIZE,
                                        2 * TILE_SIZE);
                        threadgroup_barrier(mem_flags::mem_threadgroup);
                        const uint col = j * 2 + col_bits;
                        const uint out_col = n_blk * SIMD_TILES_N * 2 + col;
                        const uint out_idx = out_col * inout_stride + chunk_idx;
                        const float g = tanh_fast(simd_res_buf[sid][rb_idx + 0]);
                        const float i = sigmoid(simd_res_buf[sid][rb_idx + 1]);
                        const float f = sigmoid(simd_res_buf[sid][rb_idx + 2]);
                        const float o = sigmoid(simd_res_buf[sid][rb_idx + 3]);
                        const float state = f * state_buf[out_idx] + i * g;
                        const float h = o * tanh_fast(state);
                        state_buf[out_idx] = state;
                        simd_out_buf[sid][row * TILE_SIZE + col] = h;
                    }
                    threadgroup_barrier(mem_flags::mem_threadgroup);
                    simdgroup_load(mm.A, simd_out_buf[sid], TILE_SIZE);
                    simdgroup_store(mm.A,
                                    (n_blk < n_blks - int(simdgroups)) ? temp_result_buf : out,
                                    inout_stride, ulong2(out_chunk_base, n_blk * TILE_SIZE));
                }
            }
        }
        threadgroup_barrier(mem_flags::mem_device);
        for (int m_blk = gid; m_blk < m_blks; m_blk += threadgroups) {
            for (int n_blk = sid; n_blk < n_blks - int(simdgroups); n_blk += simdgroups) {
                for (int i = 0; i < SIMD_TILES_M; ++i) {
                    uint out_chunk_base = (m_blk * SIMD_TILES_M + i) * TILE_SIZE;
                    simdgroup_load(mm.A, temp_result_buf, inout_stride,
                                   ulong2(out_chunk_base, n_blk * TILE_SIZE));
                    simdgroup_store(mm.A, out, inout_stride,
                                    ulong2(out_chunk_base, n_blk * TILE_SIZE));
                }
            }
        }
    }
}

struct LinearArgs {
    int in_batch_tiles;
    int in_batch_tile_offset;
    int out_batch_tiles;
    int chunk_size;
};

kernel void linear(
        device const LinearArgs* const args,
        device ftype* const in_buf,
        device const ftype* const weights_buf,
        device void* const out_buf,
        // The size of this buffer is set via MTL::ComputeCommandEncoder.
        // It depends on the SIMD group count.
        threadgroup ftype (* const simd_out_buf)[TILE_SIZE * TILE_SIZE],
        KERNEL_INDEX_INPUTS)
{
    const int chunk_size = args->chunk_size;
    const int in_batch_tiles = args->in_batch_tiles;
    const int in_batch_tile_offset = args->in_batch_tile_offset;
    const int out_batch_tiles = args->out_batch_tiles;
    const int m_blks = out_batch_tiles / SIMD_TILES_M;
    const int n_blks = kLinearInnerDim / (TILE_SIZE * SIMD_TILES_N);
    const int k_tiles = kLinearContractDim / TILE_SIZE;
    const int in_stride = in_batch_tiles * TILE_SIZE;
    const int w_stride = kLinearInnerDim;
    const int out_stride = kLinearInnerDim;
    device const ftype* const bias = weights_buf + kLinearContractDim * w_stride;
    MatMul<SIMD_TILES_M, SIMD_TILES_N, 1> mm;

    for (int ts = gid; ts < chunk_size; ts += threadgroups) {
        auto in = in_buf + in_batch_tile_offset * TILE_SIZE +
                  (ts + 1) * in_batch_tiles * TILE_SIZE * kLinearContractDim;
        const auto out_buf_offset = ts * kLinearInnerDim * out_batch_tiles * TILE_SIZE;
        device auto* const out_int8 = (device int8_t*)out_buf + out_buf_offset;
        device auto* const out_ftype = (device ftype*)out_buf + out_buf_offset;

        for (int m_blk = 0; m_blk < m_blks; ++m_blk) {
            for (int n_blk = sid; n_blk < n_blks; n_blk += simdgroups) {
                mm.mm_bias(0, k_tiles, in, in_stride, m_blk * SIMD_TILES_M * TILE_SIZE, 0,
                           weights_buf, w_stride, n_blk * SIMD_TILES_N * TILE_SIZE, 0, bias);
                for (int i = 0; i < SIMD_TILES_M; ++i) {
                    for (int j = 0; j < SIMD_TILES_N; ++j) {
                        // Store this 8x8 tile to threadgroup memory as ftype.
                        simdgroup_store(mm.accum[i][j], simd_out_buf[sid], TILE_SIZE);
                        
                        const uint tile_i = (m_blk * SIMD_TILES_M + i) * TILE_SIZE;
                        const uint tile_j = (n_blk * SIMD_TILES_N + j) * TILE_SIZE;

                        // Apply tanh activation or clamping, scaling, and type conversion.
                        // Store to the output buffer.
                        for (int elem = tid & 31; elem < TILE_SIZE * TILE_SIZE; elem += 32) {
                            const ftype matmul_output = simd_out_buf[sid][elem];
                            const auto with_clamp = kLinearOutputClamp ? clamp(matmul_output, ftype(-5.0f), ftype(5.0f)) : matmul_output;
                            const auto with_tanh = kLinearOutputTanh ? tanh_fast(with_clamp) : with_clamp;
                            const auto with_scale = with_tanh * kLinearOutputScale;

                            const int in_tile_i = elem / TILE_SIZE;
                            const int in_tile_j = elem % TILE_SIZE;
                            if (kLinearOutputAsByte) {
                                out_int8[(tile_i + in_tile_i) * out_stride + tile_j + in_tile_j] =
                                    static_cast<int8_t>(with_scale);
                            } else {
                                out_ftype[(tile_i + in_tile_i) * out_stride + tile_j + in_tile_j] =
                                    static_cast<ftype>(with_scale);
                            }
                        }
                    }
                }
            }
        }
    }
}
