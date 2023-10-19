# incomes, assets, debts
# high net worth oversample
# pig bank laproscope
scf_MIcombine <-
	function (results, variances, call = sys.call(), df.complete = Inf, ...) {
		m <- length(results)
		oldcall <- attr(results, "call")
		if (missing(variances)) {
			variances <- suppressWarnings(lapply(results, vcov))
			results <- lapply(results, coef)
		}
		vbar <- variances[[1]]
		cbar <- results[[1]]
		for (i in 2:m) {
			cbar <- cbar + results[[i]]
			# MODIFICATION:
			# vbar <- vbar + variances[[i]]
		}
		cbar <- cbar/m
		# MODIFICATION:
		# vbar <- vbar/m
		evar <- var(do.call("rbind", results))
		r <- (1 + 1/m) * evar/vbar
		df <- (m - 1) * (1 + 1/r)^2
		if (is.matrix(df)) df <- diag(df)
		if (is.finite(df.complete)) {
			dfobs <- ((df.complete + 1)/(df.complete + 3)) * df.complete *
			vbar/(vbar + evar)
			if (is.matrix(dfobs)) dfobs <- diag(dfobs)
			df <- 1/(1/dfobs + 1/df)
		}
		if (is.matrix(r)) r <- diag(r)
		rval <- list(coefficients = cbar, variance = vbar + evar *
		(m + 1)/m, call = c(oldcall, call), nimp = m, df = df,
		missinfo = (r + 2/(df + 3))/(r + 1))
		class(rval) <- "MIresult"
		rval
	}
library(haven)

scf_dta_import <-
	function( this_url ){
		
		this_tf <- tempfile()
		
		download.file( this_url , this_tf , mode = 'wb' )
		
		this_tbl <- read_dta( this_tf )
		
		this_df <- data.frame( this_tbl )
		
		file.remove( this_tf )
		
		names( this_df ) <- tolower( names( this_df ) )
		
		this_df
	}

scf_df <- scf_dta_import( "https://www.federalreserve.gov/econres/files/scf2022s.zip" )

ext_df <- scf_dta_import( "https://www.federalreserve.gov/econres/files/scfp2022s.zip" )

scf_rw_df <- scf_dta_import( "https://www.federalreserve.gov/econres/files/scf2022rw1s.zip" )

stopifnot( nrow( scf_df ) == nrow( scf_rw_df ) * 5 )
stopifnot( nrow( scf_df ) == nrow( ext_df ) )
stopifnot( all( sort( intersect( names( scf_df ) , names( ext_df ) ) ) == c( 'y1' , 'yy1' ) ) )
stopifnot( all( sort( intersect( names( scf_df ) , names( scf_rw_df ) ) ) == c( 'y1' , 'yy1' ) ) )
stopifnot( all( sort( intersect( names( ext_df ) , names( scf_rw_df ) ) ) == c( 'y1' , 'yy1' ) ) )
scf_rw_df[ , 'y1' ] <- NULL

scf_df[ , 'five' ] <- 5
# scf_fn <- file.path( path.expand( "~" ) , "SCF" , "this_file.rds" )
# saveRDS( scf_df , file = scf_fn , compress = FALSE )
# scf_df <- readRDS( scf_fn )
library(stringr)

s1_df <- scf_df[ str_sub( scf_df[ , 'y1' ] , -1 , -1 ) == 1 , ]
s2_df <- scf_df[ str_sub( scf_df[ , 'y1' ] , -1 , -1 ) == 2 , ]
s3_df <- scf_df[ str_sub( scf_df[ , 'y1' ] , -1 , -1 ) == 3 , ]
s4_df <- scf_df[ str_sub( scf_df[ , 'y1' ] , -1 , -1 ) == 4 , ]
s5_df <- scf_df[ str_sub( scf_df[ , 'y1' ] , -1 , -1 ) == 5 , ]
scf_imp <- list( s1_df , s2_df , s3_df , s4_df , s5_df )

scf_list <- lapply( scf_imp , merge , ext_df )

scf_rw_df[ is.na( scf_rw_df ) ] <- 0

scf_rw_df[ , paste0( 'wgt' , 1:999 ) ] <-
	scf_rw_df[ , paste0( 'wt1b' , 1:999 ) ] * scf_rw_df[ , paste0( 'mm' , 1:999 ) ]

scf_rw_df <- scf_rw_df[ , c( 'yy1' , paste0( 'wgt' , 1:999 ) ) ]
scf_list <- lapply( scf_list , function( w ) w[ order( w[ , 'yy1' ] ) , ] )

scf_rw_df <- scf_rw_df[ order( scf_rw_df[ , 'yy1' ] ) , ]
library(survey)
library(mitools)

scf_design <- 
	svrepdesign( 
		weights = ~wgt , 
		repweights = scf_rw_df[ , -1 ] , 
		data = imputationList( scf_list ) , 
		scale = 1 ,
		rscales = rep( 1 / 998 , 999 ) ,
		mse = FALSE ,
		type = "other" ,
		combined.weights = TRUE
	)
	
scf_design <- 
	update( 
		scf_design , 
		
		hhsex = factor( hhsex , levels = 1:2 , labels = c( "male" , "female" ) ) ,
		
		married = as.numeric( married == 1 ) ,
		
		edcl = 
			factor( 
				edcl , 
				levels = 1:4 ,
				labels = 
					c( 
						"less than high school" , 
						"high school or GED" , 
						"some college" , 
						"college degree" 
					) 
			)

	)
scf_MIcombine( with( scf_design , svyby( ~ five , ~ five , unwtd.count ) ) )

scf_MIcombine( with( scf_design , svyby( ~ five , ~ hhsex , unwtd.count ) ) )
scf_MIcombine( with( scf_design , svytotal( ~ five ) ) )

scf_MIcombine( with( scf_design ,
	svyby( ~ five , ~ hhsex , svytotal )
) )
scf_MIcombine( with( scf_design , svymean( ~ networth ) ) )

scf_MIcombine( with( scf_design ,
	svyby( ~ networth , ~ hhsex , svymean )
) )
scf_MIcombine( with( scf_design , svymean( ~ edcl ) ) )

scf_MIcombine( with( scf_design ,
	svyby( ~ edcl , ~ hhsex , svymean )
) )
scf_MIcombine( with( scf_design , svytotal( ~ networth ) ) )

scf_MIcombine( with( scf_design ,
	svyby( ~ networth , ~ hhsex , svytotal )
) )
scf_MIcombine( with( scf_design , svytotal( ~ edcl ) ) )

scf_MIcombine( with( scf_design ,
	svyby( ~ edcl , ~ hhsex , svytotal )
) )
scf_MIcombine( with( scf_design ,
	svyquantile(
		~ networth ,
		0.5 , se = TRUE , interval.type = 'quantile' 
) ) )

scf_MIcombine( with( scf_design ,
	svyby(
		~ networth , ~ hhsex , svyquantile ,
		0.5 , se = TRUE , interval.type = 'quantile' ,
		ci = TRUE 
) ) )
scf_MIcombine( with( scf_design ,
	svyratio( numerator = ~ income , denominator = ~ networth )
) )
sub_scf_design <- subset( scf_design , lf == 1 )
scf_MIcombine( with( sub_scf_design , svymean( ~ networth ) ) )
this_result <-
	scf_MIcombine( with( scf_design ,
		svymean( ~ networth )
	) )

coef( this_result )
SE( this_result )
confint( this_result )
cv( this_result )

grouped_result <-
	scf_MIcombine( with( scf_design ,
		svyby( ~ networth , ~ hhsex , svymean )
	) )

coef( grouped_result )
SE( grouped_result )
confint( grouped_result )
cv( grouped_result )
degf( scf_design$designs[[1]] )
scf_MIcombine( with( scf_design , svyvar( ~ networth ) ) )
# SRS without replacement
scf_MIcombine( with( scf_design ,
	svymean( ~ networth , deff = TRUE )
) )

# SRS with replacement
scf_MIcombine( with( scf_design ,
	svymean( ~ networth , deff = "replace" )
) )
# MIsvyciprop( ~ married , scf_design ,
# 	method = "likelihood" )
# MIsvyttest( networth ~ married , scf_design )
# MIsvychisq( ~ married + edcl , scf_design )
glm_result <- 
	scf_MIcombine( with( scf_design ,
		svyglm( networth ~ married + edcl )
	) )
	
summary( glm_result )
mean_net_worth <- scf_MIcombine( with( scf_design , svymean( ~ networth ) ) )

stopifnot( round( coef( mean_net_worth ) / 1000 , 2 ) == 1059.50 )
stopifnot( abs( 23.2 - round( SE( mean_net_worth ) / 1000 , 1 ) ) < 0.1 )
# compute quantile with all five implicates stacked (not the recommended technique)
fake_design <- svydesign( ~ 1 , data = ext_df[ c( 'networth' , 'wgt' ) ] , weights = ~ wgt )

median_net_worth_incorrect_errors <- svyquantile( ~ networth , fake_design , 0.5 )

stopifnot( round( coef( median_net_worth_incorrect_errors ) / 1000 , 2 ) == 192.7 )
library(convey)
scf_design$designs <- lapply( scf_design$designs , convey_prep )

scf_MIcombine( with( scf_design , svygini( ~ networth ) ) )
