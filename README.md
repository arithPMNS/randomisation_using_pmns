# randomisation_using_pmns
Contains examples of codes for randomised modular and scalar multiplications.
<br />
<br />
<br />
'modular_multiplication' directory contains three types of modular multiplication using the PMNS which parameters are given below:
 - p = 87798663188023528073030994638480388226916752551484470548433013409422826400553 (a 256-bit prime number)
 - n = 6
 - z = 7
 - rho = 2^50
 - gamma = 70543026063874968083330835012161440042437135154375017765226926371162653011399
 - E(X) = X^6 - 3 
 - M(X) = 1565749579740.X^5 + 123136832359.X^4 − 1697278502061.X^3 + 2943623922815.X^2 − 838163884223.X + 2587093525133
 
 The subdirectory 'mult_T0' contains number of cycles measurement for non-randomised modular multiplication.
 
 The subdirectory 'mult_T1' contains number of cycles measurement for randomised modular multiplication, where randomisation polynomial is generated BEFORE calling multiplication function.
 
 The subdirectory 'mult_T2' contains number of cycles measurement for randomised modular multiplication, where randomisation polynomial is generated IN multiplication function.
<br />
<br />
<br /> 
 'scalar_multiplication' directory contains five types of scalar multiplication for the elliptic curve which parameters are given below:
  - p = 2^256 - 1539 (a 256-bit prime number)
 - Curve equation : y^2 = x^3 + 5
Scalar multiplications are done using Algorithm 9 in https://www.matthieurivain.com/files/jcen11b.pdf 

The subdirectory 'coZ_mont_kP_s0' contains number of cycles measurement for non-randomised scalar multiplication.

The subdirectory 'coZ_mont_kP_s1' contains number of cycles measurement for non-randomised scalar multiplication. However, base point conversion (from binary to PMNS) is randomised.

The subdirectory 'coZ_mont_kP_s2' contains number of cycles measurement for randomised scalar multiplication. Base point conversion is also randomised.
Here, the same randomisation polynomial is used for the conversion process and all the scalar multiplication.

The subdirectory 'coZ_mont_kP_s3' contains number of cycles measurement for randomised scalar multiplication. Base point conversion is also randomised.
Here, the same randomisation polynomial is used for each bit of the scalar, during the multiplication scalar.

The subdirectory 'coZ_mont_kP_s4' contains number of cycles measurement for randomised scalar multiplication. Base point conversion is also randomised.
Here, a new randomisation polynomial is generated each time a modular multiplication or modular squaring is performed.
