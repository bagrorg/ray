add_library(vector vector.c)

add_library(common scene.c)
target_link_libraries(common cglm vector)

