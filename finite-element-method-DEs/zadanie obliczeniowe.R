a <- 0
b <- 2
n <- 1000
h <- 2 / (n - 1)


x_i <- function(i) {
  return((i - 1) * h)
}

e <- function(i, x) {
  if (x > x_i(i-1) && x < x_i(i)){
    return (1/h*(x - x_i(i-1)))
  }
  else if (x >= x_i(i) && x < x_i(i+1)){
    return (-1/h*(x - x_i(i+1)))
  }
  else {
    return (0)
  }
}

deriv_e <- function(i, x){
  if (x > x_i(i-1) && x < x_i(i)){
    return (1/h)
  }
  else if (x >= x_i(i) && x < x_i(i+1)){
    return (-1/h)
  }
  else {
    return (0)
  }
}


gaussian_quadrature_2p <- function(f, a, b) {
  g <- function(t) {
    f((b - a) / 2 * t + (a + b) / 2)
  }
  
  integral <- (b-a)/2*(g(-1/sqrt(3)) + g(1/sqrt(3)))
  return(integral)
}

plot_basis_functions <- function() {
  x_vals <- seq(0, 2, length.out = 200)
  plot(NULL, xlim = c(0, 2), ylim = c(0, 1), xlab = "x", ylab = "", main = paste("Basis Functions for n =", n))
  
  for (i in 1:n) {
    y_vals <- sapply(x_vals, e, i=i)
    lines(x_vals, y_vals, col = i)
  }
}
plot_basis_functions()


B <- function(i, j) {
  under_integral <-  function(x) {
    return (deriv_e(i,x)*deriv_e(j,x))
  }
  a = max(0, x_i(i-1), x_i(j-1))
  b = min(x_i(i+1), x_i(j+1))

  return (-gaussian_quadrature_2p(under_integral, a,b) + e(i, 0)*e(j,0))
}

L <- function(j) {
  return (17*e(j, 0))
}

create_A_matrix <- function() {
  A <- matrix(0, nrow = n, ncol = n)
  for (i in 1:(n - 2)) {
    A[i, i] <- B(i, i)
    A[i, i + 1] <- B(i, i + 1)
    A[i + 1, i] <- A[i, i + 1]  # Symmetry
  }
  A[n-1, n-1] = B(n-1, n-1)
  A[n, n] <- 1
  return (A)
}

create_b_vector <- function() {
  b <- sapply(1:n, L)
  b[n]<- 0
  return (b)
}


create_u_function <- function() {
  w_i <- solve(create_A_matrix(), create_b_vector())
  
  u_tilde <- function(x) {
    return (3)
  }
  
  return (function(x) {
    total <- 0
    for(i in 1:n) {
      total <- total +  w_i[i]*e(i, x)
    }
    return (total + u_tilde(x))
  })
}


plot_f <- function(f) {
  x_vals <- seq(0, 2, by = 0.01)
  y_vals <- sapply(x_vals, f)
  
  # Create the plot
  plot(x_vals, y_vals, type = 'l', col = 'blue', lwd = 2,
       xlab = 'x', ylab = 'u(x)', main = 'Graph of u(x)',
       xlim = c(0, 2), ylim = c(0,40),
       xaxs = 'i', yaxs = 'i') # 'i' for internal axis
  
  grid(nx = NULL, ny = NULL, col = "gray", lty = "dotted")
}

u = create_u_function()
plot_f(u)

