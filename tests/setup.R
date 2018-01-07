if ( .Platform$OS.type == 'windows' ) memory.limit( 256000 )

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
lodown( "scf" , scf_cat )

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
		mse = TRUE ,
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
lodown:::scf_MIcombine( with( scf_design , svyby( ~ one , ~ one , unwtd.count ) ) )

lodown:::scf_MIcombine( with( scf_design , svyby( ~ one , ~ hhsex , unwtd.count ) ) )
lodown:::scf_MIcombine( with( scf_design , svytotal( ~ one ) ) )

lodown:::scf_MIcombine( with( scf_design ,
	svyby( ~ one , ~ hhsex , svytotal )
) )
lodown:::scf_MIcombine( with( scf_design , svymean( ~ networth ) ) )

lodown:::scf_MIcombine( with( scf_design ,
	svyby( ~ networth , ~ hhsex , svymean )
) )
lodown:::scf_MIcombine( with( scf_design , svymean( ~ edcl ) ) )

lodown:::scf_MIcombine( with( scf_design ,
	svyby( ~ edcl , ~ hhsex , svymean )
) )
lodown:::scf_MIcombine( with( scf_design , svytotal( ~ networth ) ) )

lodown:::scf_MIcombine( with( scf_design ,
	svyby( ~ networth , ~ hhsex , svytotal )
) )
lodown:::scf_MIcombine( with( scf_design , svytotal( ~ edcl ) ) )

lodown:::scf_MIcombine( with( scf_design ,
	svyby( ~ edcl , ~ hhsex , svytotal )
) )
lodown:::scf_MIcombine( with( scf_design , svyquantile( ~ networth , 0.5 , se = TRUE ) ) )

lodown:::scf_MIcombine( with( scf_design ,
	svyby( 
		~ networth , ~ hhsex , svyquantile , 0.5 ,
		se = TRUE , keep.var = TRUE , ci = TRUE 
) ) )
lodown:::scf_MIcombine( with( scf_design ,
	svyratio( numerator = ~ income , denominator = ~ networth )
) )
sub_scf_design <- subset( scf_design , lf == 1 )
lodown:::scf_MIcombine( with( sub_scf_design , svymean( ~ networth ) ) )
this_result <-
	lodown:::scf_MIcombine( with( scf_design ,
		svymean( ~ networth )
	) )

coef( this_result )
SE( this_result )
confint( this_result )
cv( this_result )

grouped_result <-
	lodown:::scf_MIcombine( with( scf_design ,
		svyby( ~ networth , ~ hhsex , svymean )
	) )

coef( grouped_result )
SE( grouped_result )
confint( grouped_result )
cv( grouped_result )
degf( scf_design$designs[[1]] )
lodown:::scf_MIcombine( with( scf_design , svyvar( ~ networth ) ) )
# SRS without replacement
lodown:::scf_MIcombine( with( scf_design ,
	svymean( ~ networth , deff = TRUE )
) )

# SRS with replacement
lodown:::scf_MIcombine( with( scf_design ,
	svymean( ~ networth , deff = "replace" )
) )
lodown:::MIsvyciprop( ~ married , scf_design ,
	method = "likelihood" )
lodown:::MIsvyttest( networth ~ married , scf_design )
lodown:::MIsvychisq( ~ married + edcl , scf_design )
glm_result <- 
	lodown:::scf_MIcombine( with( scf_design ,
		svyglm( networth ~ married + edcl )
	) )
	
summary( glm_result )
library(convey)
scf_design$designs <- lapply( scf_design$designs , convey_prep )

lodown:::scf_MIcombine( with( scf_design , svygini( ~ networth ) ) )

