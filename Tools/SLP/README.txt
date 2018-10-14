This package includes the MATLAB implementation of SLP, which can do label propagation on large-scale graphs. ATTENTION that this version is for binary tasks. We'll release mutli-class version in the future.

Make sure you've download both code.zip and data.zip and unzip them to the same folder. The workspace folder should look like:
+--workspace folder
|---+--data                Folder
    |---covtype.mat     Data file
|---+--data_ssl            Folder
    |---covtype.mat     Data file
|---+--utils               Folder
    |---Accuracy.m      
|---demo.m              Show how my codes work
|---install.m           Compilation script
|---README.txt          This File :)
|---SLP.m               SLP function
|---find2cells.cpp      Source file
|---iteration.cpp       Source file

RUN:
This package has been test with MATLAB 2015b, the mex builds with 'Microsoft Visual C++ 2015 Professional'. You will need to configure you mex first.
After configuration, run the install.m first to compile *.cpp and then run demo.m to see how my codes work.

The main function of this package is in SLP.m file. For your own data, you should provide an affinity matrix `W` (For efficiency, the matrix is organized in the form of `target x source' ) and a label vector `y` (1 for positive and -1 for negative).

For any problem with my codes, feel free to drop me a message via liangdm@lamda.nju.edu.cn. Also, I hope you to cite my IJCAI'18 paper in your publications.

De-Ming Liang
April 30, 2018