##
## Johansen Procedure
##
ca.jo <- function(x, type=c("eigen", "trace"), constant=FALSE, K=2, spec=c("longrun", "transitory"), season=NULL, dumvar=NULL, ctable=c("A1", "A2", "A3"))
{
    x <- as.matrix(x)
    colnames(x) <- make.names(colnames(x))
    type <- match.arg(type)
    spec <- match.arg(spec)
    ctable <- match.arg(ctable)
    K <- as.integer(K)
    if (K < 2) {
        stop("\nK must be at least K=2.\n")
    }
    P <- ncol(x)
    arrsel <- P
    N <- nrow(x)
    if (!is.null(season)) {
        s <- season - 1
    } else {
        s <- 0
    }
    if (N * P < P + s * P + K * P^2 + P * (P + 1)/2) 
        stop("\nInsufficient degrees of freedom.\n")
    if (P > 5) 
        warning("\nToo many variables, critical values cannot be computed.\n")
    if (!(is.null(season))) {
        dum <- (diag(season) - 1/season)[, -season]
        dums <- dum
        while (nrow(dums) < N) {
            dums <- rbind(dums, dum)
        }
        dums <- dums[1:N, ]
        if (NA %in% x) {
            idx.NA <- 1:N
            ind <- as.logical(sapply(idx.NA, function(z) sum(is.na(x[z, 
                ]) * 1)))
            ind2 <- ind * (1:N)
            dums <- dums[-ind2, ]
        }
        colnames(dums) <- paste("sd", 1:ncol(dums), sep="")
    }
    if (!(is.null(dumvar))) {
        dumvar <- as.matrix(dumvar)
        colnames(dumvar) <- make.names(colnames(dumvar))
        if(is.null(colnames(dumvar))){
          dumcols <- ncol(dumvar)
          colnames(dumvar) <- paste("exo", 1:dumcols, sep = "")
          warning("\nNo column names in 'dumvar', using prefix 'exo' instead.\n")
        }          
        if (!(nrow(dumvar) == nrow(x))) {
            stop("\nUnequal row length between dummy variables and x matrix.\n")
        }
        if (NA %in% x) {
            idx.NA <- 1:N
            ind <- as.logical(sapply(idx.NA, function(z) sum(is.na(x[z, 
                ]) * 1)))
            ind2 <- ind * (1:N)
            dumvar <- dumvar[-ind2, ]
        }
    }
    x <- na.omit(x)
    N <- nrow(x)
    Z <- embed(diff(x), K)
    Z0 <- Z[, 1:P]
    a1 <- array(c(2.816, 12.099, 18.697, 24.712, 30.774, 3.962, 
        14.036, 20.778, 27.169, 33.178, 6.936, 17.936, 25.521, 
        31.943, 38.341, 2.816, 13.338, 26.791, 43.964, 65.063, 
        3.962, 15.197, 29.509, 47.181, 68.905, 6.936, 19.31, 
        35.397, 53.792, 76.955), c(5, 3, 2))
    a2 <- array(c(6.691, 12.783, 18.959, 24.917, 30.818, 8.083, 
        14.595, 21.279, 27.341, 33.262, 11.576, 18.782, 26.154, 
        32.616, 38.858, 6.691, 15.583, 28.436, 45.248, 65.956, 
        8.083, 17.844, 31.256, 48.419, 69.977, 11.576, 21.962, 
        37.291, 55.551, 77.911), c(5, 3, 2))
    a3 <- array(c(7.563, 13.781, 19.796, 25.611, 31.592, 9.094, 
        15.752, 21.894, 28.167, 34.397, 12.74, 19.834, 26.409, 
        33.121, 39.672, 7.563, 17.957, 32.093, 49.925, 71.472, 
        9.094, 20.168, 35.068, 53.347, 75.328, 12.741, 24.988, 
        40.198, 60.054, 82.969), c(5, 3, 2))
    if (ctable == "A1") {
        cvals <- a1
    } else if (ctable == "A2") {
        cvals <- a2
    } else if (ctable == "A3") {
        cvals <- a3
    }
    if (constant) {
        if (spec == "longrun") {
            ZK <- cbind(x[-c((N - K + 1):N), ], 1)
            Lnotation <- K
          } else if (spec == "transitory") {
            ZK <- cbind(x[-N, ], 1)[K:(N - 1), ]
            Lnotation <- 1
          }
        colnames(ZK) <- c(paste(colnames(x), ".l", Lnotation, sep=""), "constant")
        Z1 <- Z[, -c(1:P)]
        temp1 <- NULL
        for(i in 1:(K-1)){
          temp <- paste(colnames(x), ".dl", i, sep="")
          temp1 <- c(temp1, temp)
        }
        colnames(Z1) <- temp1
        P <- P + 1
        idx <- 0:(P - 2)
        model <- "without linear trend and constant in cointegration"
    } else {
        Z1 <- Z[, -c(1:P)]
        Z1 <- cbind(1, Z1)
        temp1 <- NULL
        for(i in 1:(K-1)){
          temp <- paste(colnames(x), ".dl", i, sep="")
          temp1 <- c(temp1, temp)
        }
        temp1 <- c("constant", temp1)
        colnames(Z1) <- temp1
        if (spec == "longrun") {
            ZK <- x[-c((N - K + 1):N), ]
            Lnotation <- K
          }
        else if (spec == "transitory") {
            ZK <- x[-N, ][K:(N - 1), ]
            Lnotation <- 1
           }
        colnames(ZK) <- paste(colnames(x), ".l", Lnotation, sep="")
        idx <- 0:(P - 1)
        model <- "with linear trend"
      }
    N <- nrow(Z0)
    if (!(is.null(season))) {
      if(constant) {
        Z1 <- cbind(dums[-(1:K), ], Z1)
      } else {
        Z1 <- cbind(Z1[, 1], dums[-(1:K), ], Z1[, -1])
        colnames(Z1) <- c("constant", colnames(Z1)[-1])
      }
    }
    if (!(is.null(dumvar))) {
      tmp <- colnames(Z1)
      if(constant){
        Z1 <- cbind(dumvar[-(1:K), ], Z1)
        colnames(Z1) <- c(colnames(dumvar), tmp)
      } else {
        Z1 <- cbind(Z1[, 1], dumvar[-(1:K), ], Z1[, -1])
        colnames(Z1) <- c("constant", colnames(dumvar), tmp[-1])
      }
    }
    M00 <- crossprod(Z0)/N
    M11 <- crossprod(Z1)/N
    MKK <- crossprod(ZK)/N
    M01 <- crossprod(Z0, Z1)/N
    M0K <- crossprod(Z0, ZK)/N
    MK0 <- crossprod(ZK, Z0)/N
    M10 <- crossprod(Z1, Z0)/N
    M1K <- crossprod(Z1, ZK)/N
    MK1 <- crossprod(ZK, Z1)/N
    M11inv <- solve(M11)
    R0 <- Z0 - t(M01 %*% M11inv %*% t(Z1))
    RK <- ZK - t(MK1 %*% M11inv %*% t(Z1))
    S00 <- M00 - M01 %*% M11inv %*% M10
    S0K <- M0K - M01 %*% M11inv %*% M1K
    SK0 <- MK0 - MK1 %*% M11inv %*% M10
    SKK <- MKK - MK1 %*% M11inv %*% M1K
    Ctemp <- chol(SKK, pivot = TRUE)
    pivot <- attr(Ctemp, "pivot")
    oo <- order(pivot)
    C <- t(Ctemp[, oo])
    Cinv <- solve(C)
    S00inv <- solve(S00)
    valeigen <- eigen(Cinv %*% SK0 %*% S00inv %*% S0K %*% t(Cinv))
    lambda <- valeigen$values
    e <- valeigen$vector
    V <- t(Cinv) %*% e
    Vorg <- V
    V <- sapply(1:P, function(x) V[, x]/V[1, x])
    W <- S0K %*% V %*% solve(t(V) %*% SKK %*% V)
    PI <- S0K %*% solve(SKK)
    DELTA <- S00 - S0K %*% V %*% solve(t(V) %*% SKK %*% V) %*% 
        t(V) %*% SK0
    GAMMA <- M01 %*% M11inv - PI %*% MK1 %*% M11inv
    if (type == "trace") {
        type <- "trace statistic"
        teststat <- as.matrix(rev(sapply(idx, function(x) -N * 
            sum(log(1 - lambda[(x + 1):P])))))
        colnames(teststat) <- "trace"
        if (arrsel > 5) {
            cval <- NULL
        }
        else {
            cval <- round(cvals[1:arrsel, , 2], 2)
            colnames(cval) <- c("10pct", "5pct", "1pct")
            rownames(cval) <- c(paste("r <= ", (arrsel - 1):1, 
                " |", sep = ""), "r = 0  |")
        }
    }
    else if (type == "eigen") {
        type <- "maximal eigenvalue statistic (lambda max)"
        teststat <- as.matrix(rev(sapply(idx, function(x) -N * 
            log(1 - lambda[x + 1]))))
        colnames(teststat) <- "lambda max."
        if (arrsel > 5) {
            cval <- NULL
        }
        else {
            cval <- round(cvals[1:arrsel, , 1], 2)
            colnames(cval) <- c("10pct", "5pct", "1pct")
            rownames(cval) <- c(paste("r <= ", (arrsel - 1):1, 
                " |", sep = ""), "r = 0  |")
        }
    }
  
    colnames(V) <- colnames(ZK)
    rownames(V) <- colnames(ZK) 
    rownames(W) <- paste(colnames(x), ".d", sep = "")
    colnames(W) <- colnames(ZK)
    colnames(Vorg) <- colnames(V)
    rownames(Vorg) <- rownames(V)
    rownames(PI) <- rownames(W)
    colnames(PI) <- colnames(W)
    colnames(Z0) <- paste(colnames(x), ".d", sep="")
    colnames(R0) <- paste("R0", colnames(Z0), sep=".")
    colnames(RK) <- paste("RK", colnames(ZK), sep=".")
    rownames(GAMMA) <- rownames(W)
   
    new("ca.jo", x = x, Z0 = Z0, Z1 = Z1, ZK = ZK, type = type, model = model, const = constant, lag = K, P = arrsel, season = season, dumvar = dumvar, cval = cval, teststat = as.vector(teststat), lambda = lambda, Vorg = Vorg, V = V, W = W, PI = PI, DELTA = DELTA, GAMMA = GAMMA, R0 = R0, RK = RK, bp = NA, test.name = "Johansen-Procedure", spec = spec, call = match.call())  
}
