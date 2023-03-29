# Функция задания внешнего излучения --------------------------------------

# A exp(-i * k * x + phi0) = A * cos(k * x + phi0) + iA sin(k * x + phi0)
f_out <- function(collocs_matrix, A=1.0, k=1.0, orient=c(1.0, 0.0), phi0=0.0) {
    return(A * exp(-1i * k * (collocs_matrix %*% orient) + phi0))
}

# Функция рефракции на двумерной области ----------------------------------

nu_refr_circle <- function(collocs_matrix, n_refr = 1.0 + 0.0i, radius = 1.0, center = c(0, 0)) {
    return((sqrt(rowSums((collocs_matrix - center)**2)) <= radius) * (n_refr - 1.0) + 1.0)
}

nu_refr_rectangle <- function(collocs_matrix, n_refr = 1.0 + 0.0i, 
                              h_low = -0.5, h_top = 0.5, w_low = -0.5, w_top = 0.5) {
    return(
        (collocs_matrix[, 1] <= w_top) * 
            (collocs_matrix[, 1] >= w_low) * 
            (collocs_matrix[, 2] <= h_top) *
            (collocs_matrix[, 2] >= h_low) * (n_refr - 1.0) + 1.0
        )
}

# Функция двумерного ядра оператора ---------------------------------------


Hankel <- function(x, nu = 0, deg = 1) {
    return(base::besselJ(x, nu) - (1i * (-1)^((deg - 1) == 1)) * base::besselY(x, nu))
}


Green_func_2d <- function(x, y, k = 1.0) {
    dist <- sqrt(sum((x - y) ** 2))
    return(-1i/4 * Hankel(k * dist, 0, 2))
}


Green_func_3d <- function(x, y, k = 1.0) {
    dist <- sqrt(sum((x - y) ** 2))
    return(exp(1i * k * dist) / (4 * pi * dist))
}

# Функция заполнения матрицы оператора и внешнего вектора ------------------

compute_problem <- function(
        collocs_grid, 
        kernel, 
        volumes = 1.0, 
        refrs = 1.0 + 0.0i, 
        k = 1.0,
        orientation = c(1.0, 0.0), 
        phi0 = 0.0, 
        A = 1.0
) {
    matrix <- diag(0 + 0i, nrow = nrow(collocs_grid))
    G_matrix <- diag(0 + 0i, nrow = nrow(collocs_grid))
    free_out <- f_out(collocs_grid, A, k, orientation, phi0)
    for (equation in 1:(nrow(matrix) - 1)) {
        for (sum_element in (equation + 1):ncol(matrix)) {
            matrix[equation, sum_element] <- 
                kernel(
                    collocs_grid[equation, ], 
                    collocs_grid[sum_element, ],
                    k
                )  * volumes
            matrix[sum_element, equation] <- matrix[equation, sum_element]
            G_matrix[equation, sum_element] <- matrix[equation, sum_element] * (refrs[sum_element] - 1)
            G_matrix[sum_element, equation] <- matrix[equation, sum_element] * (refrs[equation] - 1)
        }
    }
    return(list(A_matrix = G_matrix * (-k**2) + diag(1, nrow = nrow(G_matrix)), vec_b = - matrix %*% free_out))
}

# Постановка области решения ----------------------------------------------

center <- c(0.0, 0.0)
len <- c(3.0, 3.0)
n_x <- n_y <- 80

x_grid <- seq(center[1] - len[1]/2, center[1] + len[1]/2, length.out = n_x + 1)
y_grid <- seq(center[2] - len[2]/2, center[2] + len[2]/2, length.out = n_y + 1)

h_x <- len[1]/n_x
h_y <- len[2]/n_y

collocs_x <- seq(center[1] - len[1]/2 + h_x/2, center[1] + len[1]/2, h_x)
collocs_y <- seq(center[2] - len[2]/2 + h_y/2, center[2] + len[2]/2, h_y)

collocs_grid <- as.matrix(expand.grid(collocs_x, collocs_y))

print(collocs_grid)

ds <- h_x * h_y

# Постановка для внешнего источника излучения -----------------------------

c <- 330            # Скорость звука в свободной среде (воздух, метр в секунду)
omega <- 2000        # Частота женского голоса в герцах
k <- omega / c
orientation <- c(0.0, 1.0)
orientation <- orientation / sqrt(sum(orientation**2))
phi0 <- 0
Amplitude <- 1.0

# refr_field <- nu_refr_circle(
#     collocs_matrix = collocs_grid, 
#     n_refr = 2.25 + 1.25i, 
#     radius = 0.75, 
#     center = c(0.0, 0.0)
# ) + 1.0

refr_field <- nu_refr_rectangle(
    collocs_matrix = collocs_grid, 
    n_refr = 2.25 + 1.25i, 
    h_low = -0.7, h_top = 0.1, w_low = -0.5, w_top = 0.5
)

# Постановка задачи -------------------------------------------------------

problem_list <- compute_problem(
    collocs_grid = collocs_grid, 
    kernel = Green_func_2d, 
    volumes = ds, 
    refrs = refr_field, 
    k = k, 
    phi0 = phi0, 
    orientation = orientation, 
    A = Amplitude
)

A_matrix <- problem_list$A_matrix
vec_b <- problem_list$vec_b

# Решение СЛАУ ------------------------------------------------------------

result <- solve(A_matrix, vec_b)

resid <- A_matrix %*% result - vec_b


image(x = collocs_x, y = collocs_y, z = matrix(Im(result), nrow = length(collocs_x), byrow=T))

filled.contour(x = collocs_x, y = collocs_y, z = matrix(Im(result), nrow = length(collocs_x), byrow=T), color.palette = heat.colors)
filled.contour(x = collocs_x, y = collocs_y, z = matrix(Re(result), nrow = length(collocs_x), byrow=T), color.palette = heat.colors)