# version for Python: "^" replaced with "**"

#=========================================================================#
# functions to generate expressions for:                                  #
# 	EDGEWORTH EXPANSIONS,                                                 #
#   UNBIASED MOMENT ESTIMATES (including products and powers of moments)  #
# E(Xbar^k1*X2bar^k2*X3bar^k3*...)                                        #
# powvect: c(k1, k2, k3, ...), k's can be 0 if needed,                    #
#          e.g. E(Xbar^k1*X3bar^k3) -> c(k1, 0, k3)                       #
#=========================================================================#

# sorting the groups in a specific way; 
# this function creates a next grouping from a given one
# from 4 3 3 2 2 1  get 4 3 3 3 1 1
one_forward <- function(vect) {
	m <- sum(vect)
	if (length(vect) == 1) return(vect)
	for (i in (length(vect) - 1):1) {
		if (i == 1) {
			vect[i] <- vect[i] + 1
			return(c(vect[1:i], rep(1, m - sum(vect[1:i]))))
		}	
		else if (vect[i - 1] > vect[i]) {	
			vect[i] <- vect[i] + 1
			return(c(vect[1:i], rep(1, m - sum(vect[1:i]))))
		}
	}
}

# creates a list of all possible groupings using our sorting rule
groups <- function(nelem) { # number of elements, sum(powvect)
	l <- list(vect <- rep(1, nelem))
	while (length(vect) > 1) {
		vect <- one_forward(vect)
		l <- c(l, list(vect))
	}
	return(l)
}	

# all the combinations for letters where "a" is Xbar, "b" is X2bar and so on
# input is a vector of powXbar (# of a's), powX2bar (# of b's), powX3bar, ...
perm_categories <- function(powvect) {
	K <- length(powvect)
	nall <- sum(powvect)
	prev_vect <- paste(rep(0, nall), collapse = "")
	spaces_left <- c(nall - c(0, cumsum(powvect)[-K]))
	for (k in 1:K) {
		ind_next <- combn(1:spaces_left[k], powvect[k])
		prev_vect <- add_letter(ind_next, prev_vect, letters[k])
	}
	return(prev_vect)
}

# helper for perm_categories()
add_letter <- function(ind_next, prev_vect, letter) {  # auxiliary for permutation
	# can do list instead, so don't need to strsplit and collapse
	next_vect <- matrix(nrow = length(prev_vect), ncol = ncol(ind_next))
	for (i in 1:length(prev_vect)) {
		for (j in 1:ncol(ind_next)) {
			old <- strsplit(prev_vect[i], split = "")[[1]]
			not_filled <- which(old == "0")  # 0's originally
			old[not_filled[ind_next[, j]]] <- letter
			next_vect[i, j] <- paste(old, collapse = "")
		}
	}
	return(as.character(next_vect))	
}

# sort letters in a string - needed to compare strings
sort_str <- function(str) {
	return(paste(sort((strsplit(str, split = "")[[1]])), collapse = ""))
}

# calculate coefficient for one grouping
# elements sorted, K is length(powvect); the biggest letter is letters[K]
# excl is a vector (grep class) to remove letters and leave only one
calculate_coef <- function(str_vect, K, excl) {  # new
	coef_numer(str_vect, K, excl)/coef_denom(str_vect)
}	

# numerator for calculate_coef()
coef_numer <- function(str_vect, K, excl) { 
	coefs <- matrix(nrow = K, ncol = length(str_vect))
	for (k in 1:K) {
		nn <- nchar(gsub(excl[k], "", str_vect))
		cn <- c(0, cumsum(nn))[- (length(nn) + 1)]  # 0 in front, remove last element
		coefs[k, ] <- choose(sum(nn) - cn, nn)
	}
	return(prod(coefs))
}

# denominator for calculate_coef()
# permutations of repeated groups
# str_vect has to be sorted - strings and groups
coef_denom <- function(str_vect) {        
	if (length(str_vect) == 1) return(1)
	repeats <- duplicated(str_vect)
	k <- 1
	multiples <- NULL
	for (i in 2:(length(str_vect))) {
		if (repeats[i]) k <- k + 1
		else {
			multiples <- c(multiples, k)
			k <- 1
		}
	}
	multiples <- c(multiples, k)
	return(prod(factorial(multiples)))		
}

# table of all the groups and their coefficients (number of occurences)
all_groups <- function(powvect) { 
	K <- length(powvect)
	m <- sum(powvect) 
	groupings <- groups(m)                            # list 
	perms     <- perm_categories(powvect)             # vector of strings
	excl <- paste("[**", letters[1:K], "]", sep = "")  # letters to exclude

	mat <- matrix("", ncol = m + 2) 
	for (g in 1:length(groupings)) {
		grp <- groupings[[g]]
		for (p in 1:length(perms)) {
			permgroup <- substring(perms[p],          # break the string             
			             c(1, cumsum(grp[-length(grp)]) + 1), cumsum(grp))
			if (any(permgroup == "a")) next           # check order
			newrow <- sort(sapply(permgroup, sort_str))  
			                          # sorted strings and sorted vector
			mat <- rbind(mat, c(newrow, rep("", m - length(newrow)), 
			             length(newrow), calculate_coef(newrow, K, excl)))			                          
		}
	}
	all_empty <- apply(mat, 2, function(s) all(s == ""))
	mat <- unique(mat[, !all_empty])[-1, ]
	if (class(mat) != "matrix") 
		mat <- matrix(mat, nrow = 1)                  # if vector
	colnames(mat) <- c(paste("group", 1:(ncol(mat) - 2)), "d", "coef")  
    return(mat)
}

# takes output of all_groups()
# subscript of of mu is the order (moment)
# combines same order combinations and adds coefficients
combine_coef <- function(char_mat) {
	if (all(char_mat == "")) return(char_mat)
	J <- ncol(char_mat)
	if (all(dim(char_mat) == c(1, 3)))
		return(data.frame(mu = str_to_moment(char_mat[1, 1]), k = 1, coef = 1))
	num_mat <- t(apply(char_mat[, 1:(J-2)], 1, str_to_moment))
	mu_vect <- apply(num_mat, 1, vect_to_one)
	coef_vect <- tapply(char_mat[, "coef"], mu_vect,     # add coefficients
	                    function(x) sum(as.numeric(x)))
	d_vect    <- tapply(char_mat[, "d"], mu_vect, unique)	
	df_d <- data.frame(mu = names(d_vect), d = d_vect)
	df_coef <- data.frame(mu = names(coef_vect), coef = coef_vect)
	return(merge(df_d, df_coef, by = "mu"))
}

# helper for combine_coef()
# order of each group, eg "aacdd" has order 13 (mu13)
str_to_moment <- function(str_vect) {
	sapply(str_vect, function(str) {
		sum(match(strsplit(str, split = "")[[1]], letters))
	})
}
# helper for combine_coef()
# representation of moments (orders) in a grouping in single string
vect_to_one <- function(vect) {  
	vect <- sort(vect[vect != 0])
	return(paste(vect, collapse = ":"))
}

#-------------------------------------------------------------------------#
#              Generate an expression string for Sage/SymPy:              #
#                    SR(expr_str) or sympify(expr_str)                    #
#-------------------------------------------------------------------------#

# inputs: vector of three elements (mu, d, and coef), 
#         sample size notation (choose_from) 
one_grouping <- function(vect3elem, smpsize) {
	if (all(vect3elem == "")) return("")
	mu <- strsplit(vect3elem[1], split = ":")
	mutab <- rev(table(mu))
	muexpr <- paste("mu", names(mutab), "**", mutab, "*", 
    	             sep = "", collapse = "")
	muexpr <- substr(muexpr, start = 1, stop = nchar(muexpr) - 1) 
    
	k <- as.numeric(vect3elem[2])
	nexpr <- smpsize
	if (k > 1) {
		for (i in 1:(k-1)) {
			nexpr <- paste(nexpr, "*(", smpsize, "-", i, ")", 
			               sep = "", collapse = "")
		}          
	}		
	coef_n_mu <- paste(vect3elem[3], "*", nexpr, "*", muexpr, " + ", 
                       sep = "", collapse = "")
	return(coef_n_mu)
}

one_combination <- function(powvect, smpsize = "n") {
	if (!length(powvect)) return("1")
	if (powvect[1] == 1 & sum(powvect[-1]) == 0) return("0") 
	if (!sum(powvect)) return("1")  # all 0's
	
	res <- all_groups(powvect)
	combined_res <- combine_coef(res)
	combined_res[, 1] <- as.character(combined_res[, 1])
	vect <- apply(combined_res, 1, one_grouping, smpsize = smpsize)
	str <- paste(vect, collapse = "")
	return(paste(" (", substr(str, start = 1, stop = nchar(str) - 2), 
                 ") / ", smpsize, "**", sum(powvect), sep = ""))
}

#=========================================================================#
#-------------------------------------------------------------------------#
#        Generate assignments E(theta^m) for Edgeworth expansions         #
#-------------------------------------------------------------------------#
#=========================================================================#

# all the necessary combinations
generateComb <- function(K) {
	M <- K + 2
	comb <- NULL
	for (m in 1:M) {
		for (j in 0:m) {
			for (k in 1:K) {
				for (i in 0:floor(k/2)) {
					for (i3 in 0:(k - 2*i)) {
						for (i4 in 0:i) {
							comb <- rbind(comb, c(j + 2*i4, i3))
						}	
					}	
				}
			}
		}
	}
	ucomb <- unique(comb)
	return(ucomb[order(ucomb[, 2], ucomb[, 1]), ])
}

# K is the number of terms of Edgeworth expansion
# Studentized mean
generateETheta <- function(K, two.smp = FALSE) {
	FUN <- if (two.smp) EThetaPow2 else EThetaPow
	assign_strs <- character(K + 2)
	for (m in 1:(K + 2)) {
		assign_strs[m] <- paste("ET", m, " = ", FUN(K, m), sep = "")
	}
	cat(paste(assign_strs, collapse = "\n"))
}

# Standardized mean
generateETheta0 <- function(K) {
	assign_strs <- character(K + 2)
	for (m in 1:(K + 2)) {
		assign_strs[m] <- paste("ET", m, " = (n/mu2)**(", m, "/2) * (",
		                        one_combination(m), ")",
		                        sep = "") 
	}
	cat(paste(assign_strs, collapse = "\n"))
}

#-------------------------------------------------------------------------#
#                               ONE-SAMPLE                                # 
#-------------------------------------------------------------------------#

# Assignments for E(Xbar^k1*X2bar^k2). Example:
# comb <- generateComb(K)
# generateAssignments(comb, "nu")  
generateAssignments <- function(comb, varname) {  # comb is a matrix, ncol = 2
	comb_strs <- apply(comb, 1, one_combination)
	assign_strs <- paste(varname, comb[, 1], comb[, 2], " = ", comb_strs, sep = "")
	cat(paste(assign_strs, collapse = "\n"))
}

# Example for generating rho's for all needed combinations:
# cat(paste(mapply(function(k, l) generateRho("rho", "nu", k, l),
#                  comb[, 1], comb[, 2]), collapse = "\n")) 
generateRho <- function(var1, var2, k, l) {  # var1 = rho, var2 = nu
	coefs <- choose(l, 0:l) 
	terms <- paste(coefs, "*", var2, k, l:0, "*mu2", "**", 0:l, sep = "")
	expr <- paste(rep(c("+", "-"), length.out = l + 1), terms, collapse = " ")
	paste(paste(var1, k, l, sep = ""), "=", substr(expr, 3, nchar(expr)))
}

# expression for E(theta^m), used in generateETheta()
EThetaPow <- function(K, m) {
	paste("n**(", m, "/2)", " * A**(-", m, "/2) * (rho", m, 0, " + ", 
	       generateSum12(K, m), ") ", sep = "")
}

# functions for EThetaPow()
aCoef <- function(m, k) {
	require(MASS)
	a <- ifelse(k == 0, 1, 1/(factorial(k)*2^k) * (-1)^k * prod(m + 2*(0:(k - 1))))
	return(fractions(a))
}
generateSum12 <- function(K, m) {
	terms <- NULL
	for (k in 1:K) {
		for (i in 0:floor(k/2)) {
			terms <- c(terms, paste(aCoef(m, k - i) * (-1)^i * choose(k - i, i),
			                        "*(B/A)**", k - i, " * rho", m + 2*i, k - 2*i, 
			                        sep = ""))
		}
	}
	return(paste(terms, collapse = " + "))
}

#-------------------------------------------------------------------------#
#                               TWO-SAMPLE                                # 
#-------------------------------------------------------------------------#

# Assignments for E(Xbar^k1*X2bar^k2). Example:
# comb <- generateComb(K)
# generateAssignments2(comb, "nu_x", "x")  
# generateAssignments2(comb, "nu_y", "y")   
generateAssignments2 <- function(comb, varname, rv) {  # comb is a matrix, ncol = 2
	comb_strs <- apply(comb, 1, one_combination)
	comb_addrv <- add2(comb_strs, rv)
	assign_strs <- paste(varname, comb[, 1], comb[, 2], " = ", comb_addrv, sep = "")
	cat(paste(assign_strs, collapse = "\n"))
}

# helper function for generateAssignments2()
add2 <- function(str, rv) {
	gsub("n", paste("n_", rv, sep = ""),
	     gsub("mu", paste("mu_", rv, sep = ""), str)) 
}

# Example for generating rho's for all needed combinations:
# cat(paste(mapply(function(k, l) generateRhoTau("rho", "nu_x", "x", k, l),
#                  comb[, 1], comb[, 2]), collapse = "\n"))   
# cat(paste(mapply(function(k, l) generateRhoTau("tau", "nu_y", "y", k, l),
#                  comb[, 1], comb[, 2]), collapse = "\n"))
generateRhoTau <- function(var1, var2, rv, k, l) {  # var1 = rho, var2 = nu_x, rv = x
	coefs <- choose(l, 0:l) 
	terms <- paste(coefs, "*", var2, k, l:0, "*", "mu_", rv, 2, "**", 0:l, sep = "")
	expr <- paste(rep(c("+", "-"), length.out = l + 1), terms, collapse = " ")
	paste(paste(var1, k, l, sep = ""), "=", substr(expr, 3, nchar(expr)))
}

# expression for E(theta^m), used in generateETheta()
EThetaPow2 <- function(K, m) {  # not vectorized
	paste("n**(", m, "/2)", "*A**(-", m, "/2)", " * (", 
	      generateSum1(K, m), ") ", sep = "")
}

# functions for EThetaPow2()
generateSum45 <- function(i, j, k, m) {  # returns a string
	terms <- NULL
	for (u in 0:(k - 2*i)) {
		for (v in 0:i) {
			terms <- c(terms, paste(choose(k - 2*i, u)*choose(i, v), "*", 
			                        "B_x**", k - i - (u + v), "*",
			                        "B_y**", u + v, "*",
			             "rho", m - j + 2*(i - v), k - 2*i - u, "*",
			             "tau", j + 2*v, u, sep = ""))             
		}
	}
	return(paste(terms, collapse = " + "))
}
generateSum23 <- function(j, K, m) {
	terms <- NULL
	for (k in 1:K) {
		for (i in 0:floor(k/2)) {
			terms <- c(terms, paste(aCoef(m, k - i) * (-1)^i * choose(k - i, i),
			             "*A**", i - k, " * (", generateSum45(i, j, k, m),
			             ") ", sep = ""))
		}
	}
	return(paste(terms, collapse = " + "))
}
generateSum1 <- function(K, m) {
	terms <- NULL
	for (j in 0:m) {
		terms <- c(terms, paste((-1)^j * choose(m, j), " * (rho", m - j, 0,
		             "*", "tau", j, 0, " + ", generateSum23(j, K, m), ") ",
		             sep = ""))
	}
	return(paste(terms, collapse = " + "))
}

#=========================================================================#












	