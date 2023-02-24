set(POD5_VERSION 0.1.4)
set(POD5_DIR pod5-${POD5_VERSION}-${CMAKE_SYSTEM_NAME})
set(POD5_REPO "https://github.com/nanoporetech/pod5-file-format")
set(POD5_INCLUDE ${DORADO_3RD_PARTY}/${POD5_DIR}/include)

if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set(LIB_DIR "lib64")
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "^aarch64*|^arm*")
      if((CMAKE_CXX_COMPILER_ID MATCHES "GNU") AND (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0))
        set(POD5_VERSION 0.1.8)
        set(POD5_URL "${POD5_REPO}/releases/download/${POD5_VERSION}/lib_pod5-${POD5_VERSION}-linux-gcc7-arm64.tar.gz")
        set(LIB_DIR "lib")
      else()
        set(POD5_URL "${POD5_REPO}/releases/download/${POD5_VERSION}/lib_pod5-${POD5_VERSION}-linux-arm64.tar.gz")
      endif()
    else()
      set(POD5_URL "${POD5_REPO}/releases/download/${POD5_VERSION}/lib_pod5-${POD5_VERSION}-linux-x64.tar.gz")
    endif()
    set(POD5_LIBRARIES
      ${DORADO_3RD_PARTY}/${POD5_DIR}/${LIB_DIR}/libpod5_format.a
      ${DORADO_3RD_PARTY}/${POD5_DIR}/${LIB_DIR}/libarrow.a
      ${DORADO_3RD_PARTY}/${POD5_DIR}/${LIB_DIR}/libjemalloc_pic.a
    )
elseif(APPLE)
    set(POD5_URL "${POD5_REPO}/releases/download/${POD5_VERSION}/lib_pod5-${POD5_VERSION}-osx-11.0-arm64.tar.gz")
    set(POD5_LIBRARIES
      ${DORADO_3RD_PARTY}/${POD5_DIR}/lib/libpod5_format.a
      ${DORADO_3RD_PARTY}/${POD5_DIR}/lib/libarrow.a
    )
elseif(WIN32)
    set(POD5_URL "${POD5_REPO}/releases/download/${POD5_VERSION}/lib_pod5-${POD5_VERSION}-win-x64.tar.gz")
    set(POD5_LIBRARIES
      ${DORADO_3RD_PARTY}/${POD5_DIR}/lib/pod5_format.lib
      ${DORADO_3RD_PARTY}/${POD5_DIR}/lib/arrow_static.lib
      bcrypt.lib
    )
endif()

download_and_extract(${POD5_URL} ${POD5_DIR})
