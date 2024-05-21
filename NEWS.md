# dlim 0.3.0

* The default method for penalization is set to REML.  
* For the `dlim` function, the argument `model_type = "standard"` is deprecated and replaced with `model_type = "nonlinear"`. This way, the `model_type` argument has options that describe the type of modification.  
* For the `model_comparison` function, the argument `null = "DLM"` is deprecated and replaced with `null = "none"`. This way, the `null` argument has options that describe the type of interaction in the null model.  
* Custom plotting examples added to the vignette.  
* Flexibility in the `plot_DLF` function to allow for changing the exposure-time values on the x-axis.  

# dlim 0.2.0

* Update to `dlim` function, removing the requirement to pass `df_m` for linear modification model

# dlim 0.1.0

* Added a `NEWS.md` file to track changes to the package.
