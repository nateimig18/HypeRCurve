# HypeRCurve
Implementation of two custom RC curves that help isolate critical curve characteristics, such as max rate, center-stick sensitivity [slope], center-stick [non-]linear width, as well as a non-linearity parameter in the 2nd curve. The curve output & input are both normalized to help isolate the dependence of the max-rate.


# Normalized Sample Output from Hyper RC Curve 1 (varying center-stick slope & center-stick linear width)
![c1_varyingslope](https://cloud.githubusercontent.com/assets/3208983/24439317/15a7c9aa-1413-11e7-9d3c-5a82aec55c87.gif)

# Normalized Sample Output from Hyper RC Curve 2 (varying non-linearity parameter & non-linear exponential width)
![curve2_m0_15_varying_np](https://cloud.githubusercontent.com/assets/3208983/24439236/a66db270-1412-11e7-97e4-afd4f4ab76b0.gif)

# "But thats requires alot of computation &/or storage in LUT"
True if you were left only with the functions expf() & powf() which both implement pade approximations requiring the unspeakable floating point division. However we are wizards and wizards use floating point dark magic, to shrink are computational overhead sacrificing a degree of precision that still meets our tolerances and/or desired end-goal of achieving 'greater'. 

There are four functions with varying degree implementations of minimax polynomials for myLog2(), myExp2(), myPowf(), & myExp(), as well as a super efficient algorithm for calculating cos() & sin() on the same arguement myCoSi() which uses a Cody & Waite range reduction algorithm that prevents mantissa truncation and increased accuracy in the floating point modulus operator. Below is the error plot, note the size of the interval [-100*pi, 100*pi] and the bounded error +-2.1e-7. 

# Error of Cordic like implementation of cosine & sine 
![mycosi_error2](https://cloud.githubusercontent.com/assets/3208983/24441520/83235646-141f-11e7-8282-187abb547da1.jpg)
![mycosi_error](https://cloud.githubusercontent.com/assets/3208983/24441408/1ca450fa-141f-11e7-8068-01b720322e38.jpg)  
 
# Error of myPowf(x, p) with varying p, with all order polynomial approximations varied for both myLog2() [Y] & myExp2() [X] 
![powf_xlog2_yexp2](https://cloud.githubusercontent.com/assets/3208983/24440423/a8b7f61a-1419-11e7-9c11-973592d1ea6f.png)
Note this should only be used for the special case where the input is normalized on [0, 1], the error starts to compound exponentially as can be seen in the above plot. Any improvements requires far to much computation. 

