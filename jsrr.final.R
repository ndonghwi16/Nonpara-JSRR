#
# joint sparse reduced rank regression final version
# for submission to Journal of Applied Statistics
#
jsrr = function(response,
                basis, 
                xi = 1,
                number.gammas = 5,
                number.lambdas = 100,
                lambda.min = 1e-4,
                lambda.max = 1,
                max.iterations = 500,
                epsilon = 1e-6,
                penalty = "both",
                wavelet = FALSE,
                verbose = T)
{
   sample.size = nrow(response); number.curves = ncol(response); dimension = ncol(basis)
   inverse.zero = solve(crossprod(basis) + xi * diag(dimension))
   basis.zero = basis
   if (wavelet == T)
   {
      number.phi = attr(basis.zero, "number.phi")
   }
   else
   {
      number.phi = 1
   }
   mu = colMeans(response)
   response = centering(response, mu)
   # gamma and lambda candidates
   switch(penalty,
          both = 
             {
                gamma.candidates = (1 : number.gammas) / (number.gammas + 1)
                proximal.type = c("group", "nuclear")
                number.proximal = length(proximal.type)
             },
          nuclear = 
             {
                gamma.candidates = 0
                number.gammas = 1
                proximal.type = c("nuclear")
                number.proximal = length(proximal.type)
             },
          group = 
             {
                gamma.candidates = 1
                number.gammas = 1
                proximal.type = c("group")
                number.proximal = length(proximal.type)
             },
          stop("Available options are nuclear, group, and both.")
          )
   lambda.min = 1e-6; lambda.max = 1
   lambda.candidates = exp(seq(log(lambda.min), log(lambda.max),
                               length = number.lambdas))
   scale.lambda = sample.size * number.curves
   lambda.candidates = scale.lambda * lambda.candidates
   jsrr = messenger = proxy = list()
   bic.list = c()
   number.fits = 0
   
   for (gamma.index in 1 : number.gammas)
   {
      gamma = gamma.candidates[gamma.index]
      if (number.proximal == 1)
      {
         ratio = 1
      }
      else
      {
         ratio = c(gamma, 1 - gamma)
      }
      for (k in 1 : number.proximal)
         messenger[[k]] = proxy[[k]] = matrix(0, dimension, number.curves)
      consensus = coefficients = u = matrix(0, nrow = dimension, ncol = number.curves)
      basis = basis.zero
      basis.index = 1 : dimension
      inverse = inverse.zero

      for (lambda.index in 1 : number.lambdas)
      {
         lambda = lambda.candidates[lambda.index] * ratio
         residuals = response - basis %*% coefficients
         rss.lambda = 0.5 * mean(residuals^2) + 
            penalty(coefficients, proximal.type, lambda, number.phi)
         rss.store = NULL
         for (iteration in 1 : max.iterations)
         {
            coefficients = inverse %*% (crossprod(basis, response) + xi * (consensus - u))
            for (k in 1 : number.proximal)
               proxy[[k]] = proximal(proximal.type[k], 
                                     consensus - messenger[[k]], 
                                     lambda[k] / xi,
                                     number.phi)
            consensus = coefficients / (number.proximal + 1)
            for (k in 1 : number.proximal)
               consensus = consensus + proxy[[k]] / (number.proximal + 1)
            u = u + coefficients - consensus
            for (k in 1 : number.proximal)
               messenger[[k]] = proxy[[k]] - consensus + messenger[[k]]
            residuals = response - basis %*% coefficients
            penalties = 0
            for (k in 1 : number.proximal)
               penalties = penalties + penalty(proxy[[k]], proximal.type[k], lambda[k], number.phi)
            rss.lambda.new = 0.5 * mean(residuals^2) + penalties / (number.lambdas * number.curves)
            
            if (abs(rss.lambda.new - rss.lambda) < epsilon)
               break
            rss.lambda = rss.lambda.new
         }
         for (k in number.proximal : 1)
            active.rows = rowSums(proxy[[k]]) != 0
         for (k in 1 : number.proximal)
            rank = sum(svd(proxy[[k]])$d > 1e-15)
         number.active.rows = sum(active.rows)
         #
         # check minimality & prune nonactive rows based on proxy
         #
         if (number.active.rows == 0 || rank == 0)
         {
            break
         }
            else
         {
            coefficients = coefficients[active.rows, , drop = F]
            u = u[active.rows, , drop = F]
            consensus = consensus[active.rows, , drop = F]
            for (k in 1 : number.proximal)
            {
               proxy[[k]] = proxy[[k]][active.rows, , drop = F]
               messenger[[k]] = messenger[[k]][active.rows, , drop = F]
            }
            basis = basis[, active.rows, drop = F]
            basis.index = basis.index[active.rows]
            inverse = inverse.update(inverse, active.rows)
         }
         #
         # bic computation and store the final fit
         #
         residuals = response - basis %*% coefficients
         for (k in 1 : number.proximal)
            rank = sum(svd(proxy[[k]])$d > 1e-15)
         rss = 0.5 * mean(residuals^2)
         switch(penalty,
                both = 
                   {
                      active.dimension = number.active.rows * rank 
                   },
                group = 
                   {
                      active.dimension = number.active.rows * log(dimension)
                   },
                nuclear = 
                   {
                      active.dimension = dimension * rank
                   })
         bic = sample.size * log(dimension) * log(number.curves) * log(2 * rss) +
            log(sample.size * number.curves) * active.dimension
         number.fits = number.fits + 1
         jsrr[[number.fits]] = list(coefficients = coefficients,
                                    basis.index = basis.index, 
                                    fitted.values = basis %*% coefficients,
                                    number.active.rows = number.active.rows,
                                    rank = rank,
                                    basis = basis)
         bic.list[number.fits] = bic
         
      } # lambda.index
   } # gamma.index
   #
   # save the number of fits and the best fit
   #
   jsrr$best = jsrr[[which.min(bic.list)]]
   jsrr$bic = bic.list
   if (verbose)
   {
      cat("\n=================================\n")
      cat(" best fit: \n")
      cat("    number of active rows = ", jsrr$best$number.active.rows, "\n")
      cat("    rank = ", jsrr$best$rank, "\n")
      cat("=================================\n")
   }
   jsrr
}
# average of a list
average = function(l)
{
   Reduce(`+`, l) / length(l)
}
# centering
centering = function(y, mu)
{
   y - tcrossprod(rep(1, nrow(y)), mu)
}
# group.norm
group.norm = function(b, number.phi = 1)
{
   group.norm.b = 0
   for (m in forward(number.phi + 1, nrow(b)))
      group.norm.b = group.norm.b + sqrt(sum(b[m, ]^2))
   group.norm.b
}
# inverse.update
inverse.update = function(s, active)
{
   remove = !active
   if (sum(remove) > 0)
   {
      d = s[remove, remove]
      g = matrix(s[remove, active], sum(remove), sum(active))
      f = solve(d, g)
      inverse.active = s[active, active] - t(g) %*% f
   }
   else
   {
      inverse.active = s
   }
   inverse.active
}
# lasso.norm
lasso.norm = function(b)
{
   sum(abs(b))
}
# nuclear.norm
nuclear.norm = function(b)
{
   sum(svd(b)$d)
}
# penalty
penalty = function(b, proximal.type, lambda, number.phi)
{
   number.proximal = length(proximal.type)
   penalty = 0
   for (k in 1 : number.proximal)
   {
      switch(
         proximal.type[k],
         group =
            {
               if (lambda[k] > 0)
                  penalty = penalty + lambda[k] * group.norm(b, number.phi)
            },
         nuclear =
            {
               if (lambda[k] > 0)
                  penalty = penalty + lambda[k] * nuclear.norm(b)
            },
         lasso = 
            {
               if (lambda[k] > 0)
                  penalty = penalty + lambda[k] * lasso.norm(b, number.phi)
            },
         stop("choose group, nuclear, lasso")
      )
   }
   penalty / number.proximal
}
# proximal
proximal = function(which.proximal, v, lambda = 0, number.phi)
{
   # when lambda <= 0 proximal operator is the identity
   if (lambda <= 0)
      return(v)
   # when lambda > 0 proximal operation is carried out
   switch(
      which.proximal,
      group =
         {
            p = v
            for (j in forward(number.phi + 1, nrow(v)))
            {
               row.norm = sqrt(sum(v[j, ]^2))
               if (row.norm > lambda)
                  p[j, ] = (1 - lambda / row.norm) * v[j, ]
               else
                  p[j, ] = 0
            }
         },
      nuclear =
         {
            s = svd(v)
            d = s$d - lambda
            d[d < 0] = 0
            p = s$u %*% diag(d, nrow = length(d)) %*% t(s$v)
         },
      lasso = 
         {
            p = v
            for (j in forward(number.phi + 1, nrow(v)))
            {
               right = v[j, ] > +lambda
               p[j, right] = v[j, right] - lambda
               left = v[j, ] < -lambda
               p[j, left] = v[j, left] + lambda
            }
         },
      stop("choose group, nuclear, lasso")
   )
   p
}
# predict.wavelet
predict.wavelet = function(fit, new.predictor)
{
   number.curves = ncol(fit$coefficients)
   type = attr(fit$basis, "type")
   order = attr(fit$basis, "order")
   max.resolution = attr(fit$basis, "max.resolution")
   new.basis = wavelet(new.predictor, order, max.resolution)[, fit$basis.index]
   cbind(1, new.basis) %*% fit$coefficients
}
#
# wavelet.box version 11
#
# wavelet
wavelet = function(x, order, max.resolution)
{
   # order = 2 => min.resolution = 2
   # order = 3, 4 => min.resolution = 3
   # scale and one 
   min.resolution = ceiling(log2(2 * order - 1))
   if (min.resolution > max.resolution)
      stop("min.resolution > max.resolution")
   phi = phi.matrix(x, order, min.resolution)
   alpha = alpha.matrix(order)
   psi = psi.matrix(x, order, min.resolution, alpha)
   for (resolution in forward(min.resolution + 1, max.resolution))
      psi = cbind(psi, psi.matrix(x, order, resolution, alpha))
   basis = cbind(phi, psi)
   dimension = ncol(basis)
   number.phi = ncol(phi)
   for (j in forward(number.phi + 1, dimension))
      basis[, j] = basis[, j] / max(abs(basis[, j]))
   attr(basis, "type") = "normalized"
   attr(basis, "order") = order
   attr(basis, "dimension") = dimension
   attr(basis, "basis.index") = 1 : dimension
   attr(basis, "min.resolution") = min.resolution
   attr(basis, "max.resolution") = max.resolution
   attr(basis, "number.phi") = number.phi
   basis
}
# alpha.matrix
alpha.matrix = function(order)
{
   degree = order - 1
   two.order = 2 * order
   resolution = ceiling(log2(2 * order - 1))
   knots.order = dyadic.knots(order, resolution)
   knots.two.order = dyadic.knots(two.order, resolution + 1)
   b = matrix(0, degree, degree)
   for (l in forward(1, degree))
      for (k in forward(1, degree))
         b[l, k] = bspline(knots.order[l + order], knots.two.order, two.order, -k + two.order)
   b.inv = solve(b)
   alpha = matrix(0, degree, degree)
   ri = rep(0, degree)
   for (i in forward(-degree, -1))
   {
      for (l in forward(1, degree))
      {
         r.il = 0
         for (k in forward(0, two.order - 2 + 2 * i))
         {
            r.il = r.il - (-1)^k * cbspline(k + 1 - 2 * i, two.order) * 
               cbspline(2 * l - k, two.order)
         }
         ri[l] = r.il
      }
      alpha.i = b.inv %*% ri
      for (k in forward(-degree, -1))
         alpha[-i, -k] = alpha.i[-k]
   }
   alpha
}
# dyadic.knots
dyadic.knots = function(order, resolution)
{
   c(rep(0, order), (1 : (2^resolution - 1)) / 2^resolution, rep(1 + 1e-15, order))
}
# phi.matrix
phi.matrix = function(x, order, resolution)
{
   knots = dyadic.knots(order, resolution)
   dimension = 2^resolution + order - 1
   phi = NULL
   for (i in 1 : dimension)
      phi = cbind(phi, bspline(x, knots, order, i))
   colnames(phi) = rep("phi", ncol(phi))
   phi
}
# psi.matrix
psi.matrix = function(x, order, resolution, alpha)
{
   #
   # set up
   #
   degree = order - 1
   two.order = 2 * order
   knots.order = dyadic.knots(order, resolution)
   knots.two.order = dyadic.knots(two.order, resolution + 1)
   y = 1 - x
   #
   # interior wavelets
   #
   interior = NULL
   for (i in forward(0, 2^resolution - two.order + 1))
   {
      psi.interior = 0
      for (k in forward(0, two.order - 2))
      {
         psi.interior = psi.interior + (-1)^k * cbspline(k + 1, two.order) *
            bspline.derivative(x, knots.two.order, two.order, order, 2 * i + k + two.order)
      }
      interior = cbind(interior, psi.interior) 
   }
   #
   # boundary wavelets
   #
   zero = one = NULL
   for (i in forward(-degree, -1))
   {
      psi.zero = psi.one = 0
      for (k in forward(-degree, -1))
      {
         psi.zero = psi.zero + alpha[-i, -k] *
            bspline.derivative(x, knots.two.order, two.order, order, k + two.order)
         psi.one = psi.one + alpha[-i, -k] *
            bspline.derivative(y, knots.two.order, two.order, order, k + two.order)
      }
      for (k in forward(0, two.order - 2 + 2 * i))
      {
         psi.zero = psi.zero + (-1)^k * cbspline(k + 1 - 2 * i, two.order) *
            bspline.derivative(x, knots.two.order, two.order, order, k + two.order)
         psi.one = psi.one + (-1)^k * cbspline(k + 1 - 2 * i, two.order) *
            bspline.derivative(y, knots.two.order, two.order, order, k + two.order)
      }
      zero = cbind(zero, psi.zero)
      one = cbind(psi.one, one)
   }
   psi.matrix = cbind(zero, interior, one)
   colnames(psi.matrix) = rep("psi", ncol(psi.matrix))
   psi.matrix
}
#
# spline
#
# bspline
bspline = function(x, knots, order, which.bspline)
{
   if (order > 1)
   {
      if (knots[which.bspline + order - 1] > knots[which.bspline])
         a = (x - knots[which.bspline]) / (knots[which.bspline + order - 1] - knots[which.bspline])
      else
         a = 0
      if (knots[which.bspline + order] > knots[which.bspline + 1])
         b = (knots[which.bspline + order] - x) / (knots[which.bspline + order] - knots[which.bspline + 1])
      else
         b = 0
      return(a * bspline(x, knots, order - 1, which.bspline) +
                b * bspline(x, knots, order - 1, which.bspline + 1))
   }
   else
   {
      bspline = rep(0, length(x))
      bspline[knots[which.bspline] <= x & x < knots[which.bspline + 1]] = 1
      return(bspline)
   }
}
# bspline derivatives
bspline.derivative = function(x, knots, order, derivative, which.bspline)
{
   if (derivative > 0)
   {
      if (knots[which.bspline + order - 1] > knots[which.bspline])
         a = (order - 1) / (knots[which.bspline + order - 1] - knots[which.bspline])
      else
         a = 0
      if (knots[which.bspline + order] > knots[which.bspline + 1])
         b = (order - 1) / (knots[which.bspline + order] - knots[which.bspline + 1])
      else
         b = 0
      return(a * bspline.derivative(x, knots, order - 1, derivative - 1, which.bspline) -
                b * bspline.derivative(x, knots, order - 1, derivative - 1, which.bspline + 1))
   }
   else
      return(bspline(x, knots, order, which.bspline))
}
# cbspline
cbspline = function(x, order)
{
   knots = 0 : order
   bspline(x, knots, order, 1)
}
# backward
backward = function(from, to)
{
   if (from >= to)
      return(from : to)
   else
      return(NULL)
}
# forward
forward = function(from, to)
{
   # if (is.null(from) | is.null(to))
   #    browser()
   if (from <= to)
      return(from : to)
   else
      return(NULL)
}
# predict
predict.jsrr = function(response, new.basis, fit)
{
   active.rows = fit$basis.index
   mu.mat = matrix(colMeans(response), nrow = nrow(new.basis), ncol = ncol(response), byrow = T)
   return(new.basis[, active.rows] %*% fit$coefficients + mu.mat)
}
# mse
mse = function(true, fit)
{
   difference = true - fit
   return(mean(difference^2))
}
# mae
mae = function(true, fit)
{
   difference = abs(true - fit)
   return(mean(difference))
}
#
# log
# last updated: May 3, 2024