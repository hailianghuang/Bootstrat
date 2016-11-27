----------------------
Copyright
----------------------
Copyright (c) 2016, Hailiang Huang, Mark J. Daly and Benjamin M. Neale, Analytic and Translational Genetics Unit, Massachusetts General Hospital, Boston MA
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.

----------------------
Synopsis
----------------------
A method to correct for population stratification in testing rare variants associated with complex disorders. More details at http://biorxiv.org/content/early/2016/08/11/068999

----------------------
Prerequisites
----------------------
PLINK (version 1.9 or higher) is required for this tutorial. It can be downloaded from https://www.cog-genomics.org/plink2/

In addition, check whether the compiled binary, 'plink-bootstrat', can be executed in your LINUX system (it was statically compiled and linked on the SURFsara with GCC version 4.7.2). 

If the binary can’t be executed in your system, make appropriate changes to the Makefile and compile the PLINK program with the Bootstrat feature in the ./plink-src/ directory. Please consult the PLINK manual (http://pngu.mgh.harvard.edu/~purcell/plink/download.shtml#compilation) if you need help in compiling the code. 

If the system and the Makefile are properly configured, the source code can be compiled by simply typing 'make' in the ./plink-src/ directory. After a successful compilation, copy the combined binary to the tutorial directory: 

cp ./plink-src/plink ./plink-bootstrat

Now you are ready to move on with the tutorial.

----------------------
Usage
----------------------
An example dataset has been provided in './dat/' (taken from the PLINK tutorial). 

1. Generate the principal components
plink2 --bfile ./dat/wgas1 --pca --out wgas1

2. Calculate the distance matrix. The example uses principal components 1 and 2, and lambda of 3. 
./generateMatrix.r wgas1.eigenval wgas1.eigenvec 1,2 3

3. Use the distance matrix to conduct association tests 
./plink-bootstrat --bfile dat/wgas1 --Pperm prob_3.txt --assoc --perm --out bootstrat

4. Inspect the result
sort -k3g bootstrat.assoc.perm |head






