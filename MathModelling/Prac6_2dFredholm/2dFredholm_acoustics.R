# Постановка области решения ----------------------------------------------

center <- c(0.0, 0.0)
len <- c(2.0, 2.0)
n_x <- n_y <- 50

x_grid <- seq(center[1] - len[1]/2, center[1] + len[1]/2, length.out = n_x + 1)
y_grid <- seq(center[2] - len[2]/2, center[2] + len[2]/2, length.out = n_y + 1)

h_x <- len[1]/n_x
h_y <- len[2]/n_y

collocs_x <- seq(center[1] - len[1]/2 + h_x/2, center[1] + len[1]/2, h_x)
collocs_y <- seq(center[2] - len[2]/2 + h_y/2, center[2] + len[2]/2, h_y)

collocs_grid <- as.matrix(expand.grid(collocs_x, collocs_y))

print(collocs_grid)

ds <- h_x * h_y

# Постановка задачи -------------------------------------------------------

c <- 330            # Скорость звука в свободной среде (воздух, метр в секунду)
omega <- 180        # Частота женского голоса в герцах
k <- omega / c
orientation <- c(0.5, 1.2)
orientation <- orientation / sqrt(sum(orientation**2))


# Функция задания внешнего излучения --------------------------------------

# A exp(-i * k * x + phi0) = A * cos(k * x + phi0) + iA sin(k * x + phi0)
f_out <- function(collocs_matrix, A=1.0, k=1.0, orient=c(1.0, 0.0), phi0=0.0) {
    return(A * exp(-1i * k * (collocs_matrix %*% orient) + phi0))
}

# Функция рефракции на двумерной области ----------------------------------

nu_refr_circle <- function(collocs_matrix, n_refr = 1.0 + 0.0i, radius = 1.0, center = c(0, 0)) {
    return((sqrt(rowSums((collocs_matrix - center)**2)) <= radius) * (n_refr - 1.0) + 1.0)
}

refr_field <- nu_refr_circle(collocs_grid, n_refr = 2.25 + 1.25i, radius = 0.3, center = c(-0.2, 0.2))


# Функция двумерного ядра оператора ---------------------------------------

Hankel <- function(x, nu = 0, deg = 1) {
    return(base::besselJ(x, nu) - (1i * (-1)^((deg - 1) == 1)) * base::besselY(x, nu))
}

Green_func_2d <- function(x, y, k = 1.0) {
    dist <- sqrt(sum((x - y) ** 2))
    return(-1i/4 * Hankel(k * dist, 0, 1))
}

Green_func_2d(x = c(0.25, 0.4), y = c(-0.5, -0.6))

# Функция заполнения матрицы оператора ------------------------------------

compute_kernel <- function(collocs_grid, kernel, volumes = 1.0, refrs = 1.0 + 0.0i, k = 1.0) {
    matrix <- diag(0 + 0i, nrow = nrow(collocs_grid))
    for (equation in 1:(nrow(matrix) - 1)) {
        for (sum_element in (equation + 1):ncol(matrix)) {
            matrix[equation, sum_element] = 
                kernel(
                    collocs_grid[equation, ], 
                    collocs_grid[sum_element, ],
                    k
                )  * volumes * (refrs[sum_element] - 1)
            matrix[sum_element, equation] = matrix[equation, sum_element]
        }
    }
    return(matrix * (-k**2))
}

compute_free_vec <- function(collocs_grid, kernel, volumes = 1.0, k = 1.0, orientation = c(1.0, 0.0), phi0 = 0, A = 1) {
    matrix <- diag(0 + 0i, nrow = nrow(collocs_grid))
    free_out <- f_out(collocs_grid, A, k, orientation, phi0)
    for (equation in 1:(nrow(matrix) - 1)) {
        for (sum_element in (equation + 1):ncol(matrix)) {
            matrix[equation, sum_element] = 
                kernel(
                    collocs_grid[equation, ], 
                    collocs_grid[sum_element, ],
                    k
                ) * volumes
            matrix[sum_element, equation] = matrix[equation, sum_element]
        }
    }
    return( - matrix %*% free_out)
}

matrix_G <- compute_kernel(
    collocs_grid = collocs_grid, 
    kernel = Green_func_2d, 
    volumes = h_x * h_y, 
    refrs = refr_field, 
    k = k
)    

vec_b <- compute_free_vec(
    collocs_grid = collocs_grid, 
    kernel = Green_func_2d, 
    volumes = h_x * h_y, 
    k = k, 
    orientation = orientation, 
    phi0 = 0, 
    A = 1
)

A <- matrix_G + diag(1, nrow = nrow(matrix_G))


# Решение СЛАУ ------------------------------------------------------------

result <- solve(A, vec_b)

resid <- A %*% result - vec_b

plot(x = collocs_grid[, 1], y = collocs_grid[, 2], col = round(Re(result) * 20) %% 20)
plot(x = collocs_grid[, 1], y = collocs_grid[, 2], col = round(Im(result) * 20) %% 20)


