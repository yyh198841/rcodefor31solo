select = function(beta) if(sum(beta!=0)!=0){return(TRUE)}else{return(FALSE)}
mspe = function(x,y) mean((x - y)^2)

algorithm = function(x,y1,y2,xnew,y2new,y1new,theta1,theta2,er=1e-3,K=24){
  theta11 = theta1; theta22 = theta2; xxnew = xnew; xx = x
  p = dim(x)[2]
  kk = 1
  scope1 = scope(xx,y1)
  result1 = predict(scope1,xxnew)
  mserr1 = rep(0,2); mserr1[1] = mserr1[2] = mspe(result1,y1new)
  beta1 = matrix(0,p,K)
  for(i in 1:p){
    ll = length(scope1$beta.best[[2]][[i]])
    beta1[i,1:ll] = scope1$beta.best[[2]][[i]]
  }
  beta1[abs(beta1)<er] = 0
  l21 = rep(0,2); l21[1] = l21[2] = sqrt(sum((theta11 - beta1)^2));l21
  l = apply(beta1,1,select);l
  theta22 = theta22[l,]; theta11 = theta11[l,]; xx = xx[,l]; xxnew = xxnew[,l]
  scope2 = scope(xx,y2)
  result2 = predict(scope2,xxnew)
  mserr2 = rep(0,2); mserr2[1] = mserr2[2] = mspe(result2,y2new)
  beta2 = matrix(0,sum(l),K)
  for(i in 1:sum(l)){
    ll = length(scope2$beta.best[[2]][[i]])
    beta2[i,1:ll] = scope2$beta.best[[2]][[i]]
  }
  beta2[abs(beta2)<er] = 0
  l22 = rep(0,2); l22[1] = l22[2] = sqrt(sum((theta22 - beta2)^2))
  while(l22[kk + 1] <= l22[kk] & l21[kk+1] <= l21[kk]){
    kk = kk + 1
    l = apply(beta2,1,select);l
    theta22 = theta22[l,]; theta11 = theta11[l,]; xx = xx[,l]; xxnew = xxnew[,l]
    scope1 = scope(xx,y1)
    result1 = predict(scope1,xxnew)
    mserr1[kk + 1] = mspe(result1,y1new);mserr1
    beta1 = matrix(0,sum(l),K)
    for(i in 1:sum(l)){
      ll = length(scope1$beta.best[[2]][[i]])
      beta1[i,1:ll] = scope1$beta.best[[2]][[i]]
    }
    beta1[abs(beta1)<er] = 0
    l21[kk + 1] = sqrt(sum((theta11 - beta1)^2));l21
    l = apply(beta1,1,select);l
    theta22 = theta22[l,]; theta11 = theta11[l,]; xx = xx[,l]; xxnew = xxnew[,l]
    scope2 = scope(xx,y2)
    result2 = predict(scope2,xxnew)
    mserr2[kk + 1] = mspe(result2,y2new)
    beta2 = matrix(0,sum(l),24)
    for(i in 1:sum(l)){
      ll = length(scope2$beta.best[[2]][[i]])
      beta2[i,1:ll] = scope2$beta.best[[2]][[i]]
    }
    beta2[abs(beta2)<er] = 0
    l22[kk + 1] = sqrt(sum((theta22 - beta2)^2))
  }
  l2 = matrix(c(l21,l22),2,byrow = TRUE)
  mserr = matrix(c(mserr1,mserr2),2,byrow = TRUE)
  return(list(beta1 = beta2, beta2 = beta2, l2 = l2, mserr = mserr))
}



algorithm_simplify = function(x,y1,y2,xnew,y2new,y1new,theta1,theta2,er=1e-4,K=24){
  theta11 = theta1; theta22 = theta2; xxnew = xnew; xx = x # initial value
  l0 = p = dim(x)[2] # initial dimension
  
  scope1 = scope(xx,y1) # original estimate
  result1 = predict(scope1,xxnew); mserr1 = mspe(result1,y1new); mserr1 # original predict 
  
  beta1 = matrix(0,p,K)
  for(i in 1:p){
    ll = length(scope1$beta.best[[2]][[i]])
    beta1[i,1:ll] = scope1$beta.best[[2]][[i]]
  }
  beta1[abs(beta1)<er] = 0
  l21 = sqrt(sum((theta11 - beta1)^2)); l21
  
  l1 = apply(beta1,1,select); l1
  if(sum(l0) == 0) l1 = l0
  
  theta22 = theta22[l1,]; theta11 = theta11[l1,]; xx = xx[,l1]; xxnew = xxnew[,l1]; l0 = l1
  
  scope2 = scope(xx,y2)
  result2 = predict(scope2,xxnew); mserr2 = mspe(result2,y2new);mserr2
  
  beta2 = matrix(0,sum(l0),K)
  for(i in 1:sum(l0)){
    ll = length(scope2$beta.best[[2]][[i]])
    beta2[i,1:ll] = scope2$beta.best[[2]][[i]]
  }
  beta2[abs(beta2)<er] = 0
  l22 = sqrt(sum((theta22 - beta2)^2)); l21
  
  l1 = apply(beta2,1,select);l1
  if(sum(l0) == 0) l1 = l0
  
  theta22 = theta22[l1,]; theta11 = theta11[l1,]; xx = xx[,l1]; xxnew = xxnew[,l1]; l0 = l1
  
  scope1 = scope(xx,y1)
  result1 = predict(scope1,xxnew); mserr1[2] = mspe(result1,y1new);mserr1
  
  beta1 = matrix(0,sum(l0),K)
  for(i in 1:sum(l0)){
    ll = length(scope1$beta.best[[2]][[i]])
    beta1[i,1:ll] = scope1$beta.best[[2]][[i]]
  }
  beta1[abs(beta1)<er] = 0
  l21[2] = sqrt(sum((theta11 - beta1)^2));l21
  
  
  l2 = c(l21,l22)
  mserr = c(mserr1,mserr2)
  
  return(list(beta1 = beta2, beta2 = beta2, l2 = l2, mserr = mserr,l1 = l1))
}


example = function(){
  if(rho == 0){
    x = UniformDesignMatrix(n, p, k)
    xnew = UniformDesignMatrix(n, p, k)
  }
  if(rho !=0){
    cov_mat = (1-rho) * diag(p) + rho * matrix(1, p, p)
    x = CorrelatedDesignMatrix(n, cov_mat, k)
    xnew = CorrelatedDesignMatrix(n, cov_mat, k)
  }
  err = rnorm(n,0,sigma)
  
  y1 = y2 = y1new = y2new = rep(0,n)
  for(i in 1:n){
    xi = xinew = matrix(0,p,k)
    for(j in 1:k){
      xi[as.integer(x[i,]) == j,j] = 1
      xinew[as.integer(xnew[i,]) == j,j] = 1
    }
    y1[i] = sum(theta1*xi) + err[i]
    y2[i] = sum(theta2*xi) + err[i]
    y1new[i] = sum(theta1*xinew)
    y2new[i] = sum(theta2*xinew)
  }
  simplify = FALSE
  if(simplify) res = algorithm_simplify(x, y1, y2, xnew, y2new, y1new, theta1, theta2)
  if(!simplify) res = algorithm(x, y1, y2, xnew, y2new, y1new, theta1, theta2)
  return(res)
}


