# PartDesign for Carpentry Compiler [Project page](https://grail.cs.washington.edu/projects/carpentrycompiler/)

Our carpentry compiler converts high-level geometric designs made by users to low-level fabrication instructions that can be directly followed to manufacture parts. The compiler performs multi-objective optimization on the low-level instructions to generate Pareto-optimal candidates. 

*Note that the current implementation is developed for research. We could not gurrentee its robustness. For any problem, please email Haisen Zhao (haisen@cs.washington.edu). 

# Usage
1, Download "FreeCADLibs_12.1.4_x64_VC15" from https://github.com/FreeCAD/FreeCAD/releases/tag/0.19_pre;   
Extract the package to the local path, "D:\Develop\FreeCADLibs_12.1.4_x64_VC15".

2, Clone FreeCAD with "git clone https://github.com/chenming-wu/FreeCAD.git"  to the local path, "D:\Develop\FreeCAD";

3, Clone this repository by "git clone https://github.com/helm-compiler/carpentry-compiler.git" to the local path  "D:\Develop\Testing";  
Replace the content of "D:\Develop\FreeCAD\src\Mod\PartDesign" with "D:\Develop\Testing". Now we have only tried to compiler our framework on Windows operating system.

4, Build FreeCAD code with CMake:

    4.1 Set FREECAD_LIBPACK_DIR as the path of "FreeCADLibs_12.1.4_x64_VC15";

    4.2 Checked "BUILD_QT5";

    4.2 Configure and generate;
 
    4.3 Build to the local path, "D:\Develop\FreeCAD\build\";

5, Copy the following code into "D:\Develop\FreeCAD\src\Gui\ApplicationPy.cpp" after line #1090;
       
	if (icon.isNull())  
	{  
		QString file = QString::fromUtf8(content.c_str());  
		icon.load(file);  
		if (icon.isNull()) {  
			// ... or the name of another icon?  
			icon = BitmapFactory().pixmap(file.toUtf8());  
		}  
	}

6, Compiler FreeCAD with Visual Studio;

7, Copy some files into building folder:

     7.1. Copy the whole folder of "D:\Develop\Testing\Ext\asm3" into "D:\Develop\FreeCAD\build\Ext\freecad";  
          *"asm3" is a python-based mode for assembling;   
          https://github.com/realthunder/FreeCAD_assembly3

     7.2. Copy "carpentry_geom.dll" and "libgmp-10.dll" into "D:\Develop\FreeCAD\build\bin";  
          *"carpentry_geom.dll" includes some basic geometric functions using CGAL, et al.
         
          If you need new geometric functions, you could build the project of "carpentry-geom-lib",  
	  https://github.com/helm-compiler/carpentry-geom-lib;
	     
     7.3. Copy all files in "E:\Develop\FreeCADLibs_12.1.4_x64_VC15\bin" into "D:\Develop\FreeCAD\build\bin";
     
     7.4 Copy the folder of "E:\Develop\FreeCADLibs_12.1.4_x64_VC15\plugins\platforms" into "D:\Develop\FreeCAD\build\bin";

8, Run;


# Citation
If you use this code for your research, please cite our [paper](hhttps://grail.cs.washington.edu/projects/carpentrycompiler/files/CarpentryCompiler.pdf):

```
@article {wu_siga19,
    author = {Chenming Wu and Haisen Zhao and Chandrakana Nandi and Jeffrey I. Lipton and Zachary Tatlock and Adriana Schulz},
    title = {Carpentry Compiler},
    journal = {ACM Transactions on Graphics},
    note = {presented at SIGGRAPH Asia 2019},
    volume = {38},
    number = {6},
    pages = {Article No. 195},
    year = {2019}
}
```

# License
All rights about the program are reserved by the authors of this project. The programs can only be used for research purpose. In no event shall the author be liable to any party for direct, indirect, special, incidental, or consequential damage arising out of the use of this program.

