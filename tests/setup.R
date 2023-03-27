# incomes, assets, debts
# high net worth oversample
# pig bank laproscope

# variant of code{mitools::scf_MIcombine} that only uses the sampling variance from the first implicate instead of averaging all five

scf_scf_MIcombine <-
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

tf_s <- tempfile()

tf_p <- tempfile()

tf_rw <- tempfile()

download.file( "https://www.federalreserve.gov/econres/files/scf2019s.zip" , tf_s , mode = 'wb' )

download.file( "https://www.federalreserve.gov/econres/files/scfp2019s.zip" , tf_p , mode = 'wb' )

download.file( "https://www.federalreserve.gov/econres/files/scf2019rw1s.zip" , tf_rw , mode = 'wb' )

s_tbl <- read_dta( tf_s )

p_tbl <- read_dta( tf_p )

rw_tbl <- read_dta( tf_rw )

s_df <- data.frame( s_tbl )

p_df <- data.frame( p_tbl )

scf_rw_df <- data.frame( rw_tbl )

names( s_df ) <- tolower( names( s_df ) )

names( p_df ) <- tolower( names( p_df ) )

names( scf_rw_df ) <- tolower( names( scf_rw_df ) )

stopifnot( nrow( s_df ) == nrow( scf_rw_df ) * 5 )
stopifnot( nrow( s_df ) == nrow( p_df ) )
# scf_fn <- file.path( path.expand( "~" ) , "SCF" , "this_file.rds" )
# saveRDS( scf_df , file = scf_fn , compress = FALSE )
# scf_df <- readRDS( scf_fn )
library(survey)
library(mitools)

# confirm that the only overlapping columns
# between the three data sets are `y1`
# (the unique primary economic unit id - peu)
# and `yy1` (the five records of the peu)
stopifnot( all.equal( sort( intersect( names( s_df ) , names( p_df ) ) ) , c( 'y1' , 'yy1' ) ) )
stopifnot( all.equal( sort( intersect( names( s_df ) , names( scf_rw_df ) ) ) , c( 'y1' , 'yy1' ) ) )
stopifnot( all.equal( sort( intersect( names( p_df ) , names( scf_rw_df ) ) ) , c( 'y1' , 'yy1' ) ) )

# throw out the unique identifiers ending with `1`
# because they only match one-fifth of the records in the survey data
scf_rw_df$y1 <- NULL

# `s_df` currently contains
# five records per household -- all five of the implicates.

# add a column `one` to every record, containing just the number one
s_df$one <- 1

# add a column `five` to every record, containing just the number five
s_df$five <- 5
# note: this column should be used to calculate weighted totals.

# break `s_df` into five different data sets
# based on the final character of the column 'y1'
# which separates the five implicates
s1_df <- s_df[ substr( s_df$y1 , nchar( s_df$y1 ) , nchar( s_df$y1 ) ) == 1 , ]
s2_df <- s_df[ substr( s_df$y1 , nchar( s_df$y1 ) , nchar( s_df$y1 ) ) == 2 , ]
s3_df <- s_df[ substr( s_df$y1 , nchar( s_df$y1 ) , nchar( s_df$y1 ) ) == 3 , ]
s4_df <- s_df[ substr( s_df$y1 , nchar( s_df$y1 ) , nchar( s_df$y1 ) ) == 4 , ]
s5_df <- s_df[ substr( s_df$y1 , nchar( s_df$y1 ) , nchar( s_df$y1 ) ) == 5 , ]

# count the total number of records in `s_df`
m.rows <- nrow( s_df )

# confirm that the number of records did not change
stopifnot(
	sum( nrow( s1_df ) , nrow( s2_df ) , nrow( s3_df ) , nrow( s4_df ) , nrow( s5_df ) ) == m.rows
)

# sort all five implicates by the unique identifier
# s1_df <- s1_df[ order( s1_df$yy1 ) , ]
# s2_df <- s2_df[ order( s2_df$yy1 ) , ]
# s3_df <- s3_df[ order( s3_df$yy1 ) , ]
# s4_df <- s4_df[ order( s4_df$yy1 ) , ]
# s5_df <- s5_df[ order( s5_df$yy1 ) , ]

scf_imp <- list( s1_df , s2_df , s3_df , s4_df , s5_df )

scf_list <- lapply( scf_imp , merge , p_df )

scf_list <- lapply( scf_list , function( w ) w[ order( w[ , 'yy1' ] ) , ] )

# replace all missing values in the replicate weights table with zeroes..
scf_rw_df[ is.na( scf_rw_df ) ] <- 0

# ..then multiply the replicate weights by the multiplication factor
scf_rw_df[ , paste0( 'wgt' , 1:999 ) ] <- scf_rw_df[ , paste0( 'wt1b' , 1:999 ) ] * scf_rw_df[ , paste0( 'mm' , 1:999 ) ]

# only keep the unique identifier and the final (combined) replicate weights
scf_rw_df <- scf_rw_df[ , c( 'yy1' , paste0( 'wgt' , 1:999 ) ) ]

# sort the replicate weights data frame by the unique identifier as well
scf_rw_df <- scf_rw_df[ order( scf_rw_df$yy1 ) , ]

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
		
		hhsex = factor( hhsex , labels = c( "male" , "female" ) ) ,
		
		married = as.numeric( married == 1 ) ,
		
		edcl = 
			factor( 
				edcl , 
				labels = 
					c( 
						"less than high school" , 
						"high school or GED" , 
						"some college" , 
						"college degree" 
					) 
			)

	)
scf_MIcombine( with( scf_design , svyby( ~ one , ~ one , unwtd.count ) ) )

scf_MIcombine( with( scf_design , svyby( ~ one , ~ hhsex , unwtd.count ) ) )
scf_MIcombine( with( scf_design , svytotal( ~ one ) ) )

scf_MIcombine( with( scf_design ,
	svyby( ~ one , ~ hhsex , svytotal )
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
		0.5 , se = TRUE , method = 'constant' , interval.type = 'quantile' 
) ) )

scf_MIcombine( with( scf_design ,
	svyby(
		~ networth , ~ hhsex , svyquantile ,
		0.5 , se = TRUE , method = 'constant' , interval.type = 'quantile' ,
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

# compute mean net worth using the 2019 PUF
mean_net_worth <- scf_scf_MIcombine( with( scf_design , svymean( ~ networth ) ) )

# exactly match "Table 4" tab cell W6 of
# https://www.federalreserve.gov/econres/files/scf2019_tables_public_nominal_historical.xlsx
stopifnot( round( coef( mean_net_worth ) / 1000 , 2 ) == 746.82 )

# match mean net worth standard error within $100 from
# https://www.federalreserve.gov/publications/files/scf20.pdf#page=11
stopifnot( abs( 15.6 - round( SE( mean_net_worth ) / 1000 , 1 ) ) < 0.1 )

# compute median net worth using the 2019 PUF
median_net_worth <-
	scf_scf_MIcombine( with( scf_design ,
		svyquantile(
			~ networth ,
			0.5 , se = TRUE , method = 'constant' , interval.type = 'quantile' ,
			ci = TRUE
	) ) )

# for( this_qrule in c("math","school","shahvaish","hf1","hf2","hf3",
# 		 "hf4","hf5","hf6","hf7","hf8","hf9") ){
# 			 print( this_qrule )
# print( scf_scf_MIcombine( with( scf_design , svyquantile( ~ networth , 0.5 , interval.type = 'quantile' , qrule = this_qrule ) ) ) )
# }

# confirm the estimate and standard error match FRB calculations within one dollar
# stopifnot( round( coef( median_net_worth ) ) == 97306 )
# stopifnot( round( SE( median_net_worth ) ) == 2699 )

scf_scf_MIcombine( with( scf_design , oldsvyquantile( ~networth , 0.5 , method = 'constant' , interval.type = 'quantile' ) ) )
library(convey)
scf_design$designs <- lapply( scf_design$designs , convey_prep )

scf_MIcombine( with( scf_design , svygini( ~ networth ) ) )
