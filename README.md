# Karplus-Strong

This repository contains 5 MATLAB programs that synthesize sound using the Karplus-Strong algorithm.

KS_Basic.m is a basic implementation of the algorithm to simulate a plucked string. The parameters are fundamental frequency, loss factor, pluck strength, and duration.

KS_Tuned.m is similar to KS_Basic.m, but it makes use of an allpass filter to correct the tuning issues inherent in the basic version.

KS_Drum.m uses the modified K-S algorithm capable of producing drum-like sounds as well as other interesting sounds. It does this by introducing a blend factor parameter, which controls the probability of the sign within the algorithm.

KS_AcGuitar.m is an acoustic guitar model. A redorded pluck is used as input to the delay line and is convolved with a guitar body impulse response. A new parameter is introduced to control the pluck position along the string, which is similated by using a comb filter.

KS_ElGuitar.m is an electric guitar model. It is very similar to the acoustic guitar model, but it introduces feedback and distortion.
