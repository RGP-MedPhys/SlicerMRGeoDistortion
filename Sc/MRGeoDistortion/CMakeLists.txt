cmake_minimum_required(VERSION 2.8.9)

project(MRGeoDistortion)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "http://slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/MRGeoDistortion")
set(EXTENSION_CATEGORY "Examples")
set(EXTENSION_CONTRIBUTORS "Ryan Glen Price (Henry Ford Health System)")
set(EXTENSION_DESCRIPTION "This is an example of a simple extension")
set(EXTENSION_ICONURL "http://www.example.com/Slicer/Extensions/SlicerMRinRT.png")
set(EXTENSION_SCREENSHOTURLS "http://www.example.com/Slicer/Extensions/MRGeoDistortion/Screenshots/1.png")

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})
set(qt_module_include_directories
   ${CMAKE_CURRENT_SOURCE_DIR}/Widgets
   ${CMAKE_CURRENT_BINARY_DIR}/Widgets
   ${vtkSlicerVolumesModuleLogic_SOURCE_DIR}
   ${vtkSlicerVolumesModuleLogic_BINARY_DIR}
   ${qSlicerVolumesModuleWidgets_SOURCE_DIR}
   ${qSlicerVolumesModuleWidgets_BINARY_DIR}
   )
 set(qt_module_target_libraries
   vtkSlicerVolumesModuleLogic
   qSlicerVolumesModuleWidgets
   )

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(MeasureDistortion)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_CPACK})
