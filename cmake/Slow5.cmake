# Bring in slow5lib; it will build when something links to it.
add_subdirectory("${DORADO_3RD_PARTY_SOURCE}/slow5lib"
                 "${CMAKE_BINARY_DIR}/3rdparty/slow5lib"
                 EXCLUDE_FROM_ALL)

# Provide a stable front target that adds slow5lib src/includes to consumers.
add_library(slow5_lib INTERFACE)
target_link_libraries(slow5_lib INTERFACE slow5)

# Public headers are already exported by 'slow5', but adding them here is harmless.
# The key addition is the 'src' dir for private headers like slow5_extra.h.
target_include_directories(slow5_lib INTERFACE
    $<BUILD_INTERFACE:${DORADO_3RD_PARTY_SOURCE}/slow5lib/include>
    $<BUILD_INTERFACE:${DORADO_3RD_PARTY_SOURCE}/slow5lib/src>
)