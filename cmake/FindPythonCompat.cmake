unset(_PythonCompat_ARGS)
if(PythonCompat_FIND_VERSION)
    list(APPEND _PythonCompat_ARGS ${PythonCompat_FIND_VERSION})
endif()
if(PythonCompat_FIND_VERSION_EXACT)
    list(APPEND _PythonCompat_ARGS EXACT)
endif()
if(PythonCompat_FIND_QUIETLY)
    list(APPEND _PythonCompat_ARGS QUIET)
endif()

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.12)
    # Cmake 3.12 has FindPython.cmake, so we can simply forward to it
    unset(_PythonCompat_COMPONENTS)
    unset(_PythonCompat_OPTIONAL_COMPONENTS)
    foreach(_PythonCompat_COMPONENT ${PythonCompat_FIND_COMPONENTS})
        if(PythonCompat_FIND_REQUIRED_${_PythonCompat_COMPONENT})
            list(APPEND _PythonCompat_COMPONENTS "${_PythonCompat_COMPONENT}")
        else()
            list(APPEND _PythonCompat_OPTIONAL_COMPONENTS "${_PythonCompat_COMPONENT}")
        endif()
    endforeach()

    find_package(Python ${_PythonCompat_ARGS}
        COMPONENTS ${_PythonCompat_COMPONENTS}
        OPTIONAL_COMPONENTS ${_PythonCompat_OPTIONAL_COMPONENTS})

    set(PythonCompat_FOUND ${Python_FOUND})
    return()
endif()


if(NOT PythonCompat_FIND_COMPONENTS)
  set(PythonCompat_FIND_COMPONENTS Interpreter)
  set(PythonCompat_FIND_REQUIRED_Interpreter TRUE)
endif()

set(_PythonCompat_REQUIRED_VARS)

if(DEFINED PythonCompat_FIND_REQUIRED_Interpreter)
    if(Python_EXECUTABLE AND NOT PYTHON_EXECUTABLE)
        set(PYTHON_EXECUTABLE ${Python_EXECUTABLE} CACHE FILEPATH
            "Path to a program." FORCE)
    endif()

    find_package(PythonInterp ${_PythonCompat_ARGS})

    set(Python_Interpreter_FOUND ${PYTHONINTERP_FOUND})
    set(Python_EXECUTABLE ${PYTHON_EXECUTABLE})
    set(Python_VERSION ${PYTHON_VERSION_STRING})
    set(Python_VERSION_MAJOR ${PYTHON_VERSION_MAJOR})
    set(Python_VERSION_MINOR ${PYTHON_VERSION_MINOR})
    set(Python_VERSION_PATCH ${PYTHON_VERSION_PATCH})

    if(TARGET Python::Interpreter AND PYTHONINTERP_FOUND)
        add_executable (Python::Interpreter IMPORTED)
        set_property(TARGET Python::Interpreter PROPERTY
            IMPORTED_LOCATION "${PYTHON_EXECUTABLE}")
    endif()

    if(PythonCompat_FIND_REQUIRED_Interpreter)
        list(APPEND _PythonCompat_REQUIRED_VARS PYTHON_EXECUTABLE)
    endif()
endif()

if(DEFINED PythonCompat_FIND_REQUIRED_Development)
    find_package(PythonLibs ${_PythonCompat_ARGS})

    set(Python_Development_FOUND ${PYTHONLIBS_FOUND})
    set(Python_INCLUDE_DIRS ${PYTHON_INCLUDE_DIRS})
    set(Python_LIBRARIES ${PYTHON_LIBRARIES})
    #set(Python_LIBRARY_DIRS ${PYTHON_EXECUTABLE})
    #set(Python_RUNTIME_LIBRARY_DIRS ${PYTHON_EXECUTABLE})

    set(Python_INCLUDE_DIR ${PYTHON_INCLUDE_DIR})
    set(Python_LIBRARY_RELEASE ${PYTHON_LIBRARY_RELEASE})
    set(Python_LIBRARY_DEBUG ${PYTHON_DEBUG_LIBRARY})

    if(NOT TARGET Python::Python AND PYTHONLIBS_FOUND)
        if(PYTHON_LIBRARY MATCHES "${CMAKE_SHARED_LIBRARY_SUFFIX}$" OR
           PYTHON_DEBUG_LIBRARY MATCHES "${CMAKE_SHARED_LIBRARY_SUFFIX}$")
            set(_PythonCompat_LIBRARY_TYPE SHARED)
        else()
            set(_PythonCompat_LIBRARY_TYPE UNKNOWN)
        endif()
        add_library(Python::Python "${_PythonCompat_LIBRARY_TYPE}" IMPORTED)
        set_property(TARGET Python::Python PROPERTY
            INTERFACE_INCLUDE_DIRECTORIES "${PYTHON_INCLUDE_DIRS}")
        if(PYTHON_DEBUG_LIBRARY)
            set_property(TARGET Python::Python APPEND PROPERTY
                IMPORTED_CONFIGURATIONS RELEASE)
            set_target_properties(Python::Python PROPERTIES
                IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
                IMPORTED_LOCATION_RELEASE "${PYTHON_LIBRARY_RELEASE}")
            set_property(TARGET Python::Python APPEND PROPERTY
                IMPORTED_CONFIGURATIONS DEBUG)
            set_target_properties(Python::Python PROPERTIES
                IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "C"
                IMPORTED_LOCATION_DEBUG "${PYTHON_DEBUG_LIBRARY}")
        else()
            set_target_properties(Python::Python PROPERTIES
                IMPORTED_LINK_INTERFACE_LANGUAGES "C"
                IMPORTED_LOCATION "${PYTHON_LIBRARY}")
        endif()
    endif()

    if(PythonCompat_FIND_REQUIRED_Development)
        list(APPEND _PythonCompat_REQUIRED_VARS PYTHON_LIBRARIES
                                                PYTHON_INCLUDE_DIRS)
    endif()
endif()

include(FindPackageHandleStandardArgs)
set(Python_FIND_COMPONENTS ${PythonCompat_FIND_COMPONENTS})
find_package_handle_standard_args(Python
    REQUIRED_VARS ${_PythonCompat_REQUIRED_VARS}
    VERSION_VAR Python_VERSION
    HANDLE_COMPONENTS)
set(PythonCompat_FOUND ${Python_FOUND})
