# CMakeList.txt : CMake project for ValorCaracteristico, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

# Enable Hot Reload for MSVC compilers if supported.
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project ("ValorCaracteristico")

# Add source to this project's executable.
file(GLOB ValorCaracteristico_sources ./**.cpp)
set(ValorCaracteristico ${ValorCaracteristico_sources})
add_executable (ValorCaracteristico ${ValorCaracteristico_sources})

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET ValorCaracteristico PROPERTY CXX_STANDARD 20)
endif()

add_custom_target(results ALL
    COMMAND ${CMAKE_COMMAND} -E make_directory results)