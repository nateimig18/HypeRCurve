# HypeRCurve
Implementation of two custom RC curves that help isolate critical curve characteristics, such as max rate, center-stick sensitivity [slope], center-stick [non-]linear width, as well as a non-linearity parameter in the 2nd curve. The curve output & input are both normalized to help isolate the dependence of the max-rate.


# Normalized Sample Output from Hyper RC Curve 1 (varying center-stick slope & center-stick linear width)
![c1_varyingslope](https://cloud.githubusercontent.com/assets/3208983/24439317/15a7c9aa-1413-11e7-9d3c-5a82aec55c87.gif)

# Normalized Sample Output from Hyper RC Curve 2 (varying non-linearity parameter & non-linear exponential width)
![curve2_m0_15_varying_np](https://cloud.githubusercontent.com/assets/3208983/24439236/a66db270-1412-11e7-97e4-afd4f4ab76b0.gif)

# "But thats requires alot of computation &/or storage in LUT"
True if you were left only with the functions expf() & powf(). However we are wizards and wizards use floating point dark magic, to shrink are computational overhead sacrificing a degree of precision that still meets our tolerances and/or desired end-goal of achieving 'greater'. 
