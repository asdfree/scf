if ( .Platform$OS.type == 'windows' ) memory.limit( 256000 )

options("lodown.cachaca.savecache"=FALSE)

library(lodown)
lodown( "scf" , output_dir = file.path( getwd() ) )
library(lodown)
# examine all available SCF microdata files
scf_cat <-
	get_catalog( "scf" ,
		output_dir = file.path( getwd() ) )

# 2016 only
scf_cat <- subset( scf_cat , year == 2016 )
# download the microdata to your local computer


library(survey)
library(mitools)

scf_imp <- readRDS( file.path( getwd() , "scf 2016.rds" ) )

scf_rw <- readRDS( file.path( getwd() , "scf 2016 rw.rds" ) )

scf_design <- 
	svrepdesign( 
		weights = ~wgt , 
		repweights = scf_rw[ , -1 ] , 
		data = imputationList( scf_imp ) , 
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
		keep.var = TRUE , ci = TRUE 
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
MIsvyciprop( ~ married , scf_design ,
	method = "likelihood" )
MIsvyttest( networth ~ married , scf_design )
MIsvychisq( ~ married + edcl , scf_design )
glm_result <- 
	scf_MIcombine( with( scf_design ,
		svyglm( networth ~ married + edcl )
	) )
	
summary( glm_result )
library(convey)
scf_design$designs <- lapply( scf_design$designs , convey_prep )

scf_MIcombine( with( scf_design , svygini( ~ networth ) ) )

# compute mean net worth using the 2016 PUF
mean_net_worth <- scf_scf_MIcombine( with( scf_design , svymean( ~ networth ) ) )

# confirm the estimate and standard error match FRB calculations within one dollar
stopifnot( round( coef( mean_net_worth ) ) == 689509 )
stopifnot( round( SE( mean_net_worth ) ) == 12670 )

# compute median net worth using the 2016 PUF
median_net_worth <-
	scf_scf_MIcombine( with( scf_design ,
		svyquantile(
			~ networth ,
			0.5 , se = TRUE , method = 'constant' , interval.type = 'quantile' ,
			keep.var = TRUE , ci = TRUE
	) ) )

# confirm the estimate and standard error match FRB calculations within one dollar
stopifnot( round( coef( median_net_worth ) ) == 97306 )
stopifnot( round( SE( median_net_worth ) ) == 2699 )

