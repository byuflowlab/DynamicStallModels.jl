## Welcome To JacobChild, Jacob's Branch of Dynamic Stall Models!

A few things should be explained right off- this code is nowhere near publication ready, may contain bugs, is not optimized, and likely has lots of issues, however, it does contain some helpful changes!

This documentation will have headings as follows

- Folder
  
  - File
    
    - Line Number Explanation

Abbreviations

- DSM - Dynamic Stall Models
- f - attachment degree, the degree to which the flow is attached around the airfoil. From 0:1

#### Documentation

- **JacobsExtras**
  
  - [GettingStartedTest1.jl](JacobsExtras\GettingStartedTest1.jl)
    
    - This was to get me familiar with Dynamic Stall Models, and came from following along with the Getting Started Docs tutorial
  
  - JacobsReadMe.md
    
    - This contains notes to myself about how to manage the package and github, points 1 & 2 are outdated, as well as the plan. I learned how to use Revise.jl and develop code while in the dev folder of Julia for the DSM package. The github notes are still helpful/relevant
  
  - JacobsWeeklyUpdate.md
    
    - This was used at weekly work meetings (that I periodically attended) to update my fellows and ask advice. The images likely won't load as the path will have probably changed and some of them were temporary and not directly linked.

- **Polars**
  
  - [LarsenFig4bClStaticPlot.csv](polars\LarsenFig4bClStaticPlot.csv)
    
    - This csv is a polar created from webplot digitizing Figure 4b in Larsen's paper and extracting the Cl data.
  
  - Other files - were pre existing, or I have forgotten If I was the one to add them.

- **src**
  
  - [airfoils.jl](src\airfoils.jl)
    
    - There may be minor changes not documented, but they should have comments. Most of the function header's need to be updated
    
    - 198-200 & 256-258: This was a fix implemented to help find the max and min Cl points (used in finding separation). The max Cl is generally between 0-30degrees like we expect it, but min Cl can be outside of the -30 to 0 degree range and mess up several things downstream. Thus I search for maxcl first and then search for the mincl, by going up to that point, making sure we don't get a bad minclidx value that is after the maxclidx
    
    - 262: Around here used to be a hardcoded separation point input, I pulled it out. Call the update_airfoil function if something like that needs to be changed
    
    - 291-316: I added the update_oye_A function. This is to allow the user to use the input value found in an academic paper (say Larsen's) and convert it to the value that Hansen's function would use. This is because Hansen's tau (time delay) is hardcoded. Having this function gives the user flexibility to change their input to the proper format only if they want to. *Note: This didn't seem to make a difference in the output of f (separation degree), either implementation is wrong or it is less critical than thought.*
    
    - 656-718: The Larsen Separationpoint function (LSP). A few changes were made
      
      - 671-680: An attempt was made to make it possible to switch between Cl and Cn. It would require an additional input into the separationpoint function and thus wasn't implemented, however it would work if dsmodel was pulled in. 
      
      - 684: This line calls cl_fullsep_faber in the [Oye.jl](src\Oye.jl) file. It should be noted that there was an error in the hermite interpolation around line 249 of Oye that Weston fixed and I have input his fix as well.
      
      - 687-704: There was a numerical instability that made f  not equal to 1.0 like it should be. An if statement test was implemented to check for and correct the instability. The tolerances (tol1 and tol2, currently = .35) can be changed for the user's specific use case.
  
  - [Oye.jl](src%5COye.jl)
    
    - 19 & 74-83, 171: Version flags existed for which version (Hansen, Faber, BeddoesLeishman, Larsen) of the Oye model would be used to solve things. I added the Larsen (4) flag. Faber and Larsen have the same coefficient model, so the if statements reflect that.
    
    - 126, 193, 224, 300: I added some or all of these, I can't remember. I made sure that every hardcoded separation point equation was deleted or commented out and replaced with the multiple dispatch separationpoint function call. This allows the user to mix and match solving methods with different separationpoint equation methods.
    
    - 194-196: This is the hardcoded implementation of the hansen separationpoint function. Thus the need for the update/conversion function in airfoils.jl
    
    - 227-230: There used to be a large set of if statements to convert Tau depending on the version, this was all added and then deleted by me and the hardcoded hansen version was kept. It seems here and on line 192 are unecessary function calls to solve for cn_fs, and cn_inv (lines 192, 221, 223) these could be deleted (they are currently commented out), but it is also possible I messed up and deleted something that needs those.
    
    - 249: This is the hermite interpolation fix discussed earlier and done by Weston.

- **test**
  
  - [OyeLarsenTests.jl](test\OyeLarsenTests.jl)
    
    - This is a file meant to eventually be run as an automated github test to make sure that everything runs well and nothing was broken during other code updates. Currently it is not quite in a state to be used that way.
    
    - The code is heavily commented and should be semi easily understood from that
    
    - 21-39 & 51 - 63: This extracts the needed fdegree data and extends it a bit both directions (quite cleverly I think) to allow us to test our solved code against the outputs. It was taken from a webplot digitized polar from Fig 4c the fdyn attachment degree.
    
    - 45-50: Sets up the airfoil struct. Larsen's paper says deep stall is at 30deg, so that is hardcoded and updated. 
    
    - 64-66: Lines to plot the comparison between Larsen f and dsm f. If the test file throws an error (it currently is throwing 2) the plot will not show.
    
    - 70-88: These lines should not be used! I modified part of airfoils to output each of the cn/cl values and then I plotted and compared them all for my understanding. I then put everything back to normal. I left the code in in case I want to do it again. The resulting plot can be found [here](testing\Oye\Outputs\ClTypesPlotted.png)
    
    - 90-95: Create the airfoil struct to be tested
    
    - 105-107: Test tolerances, the first one is relative and used for the critical alpha tests, the second is absolute and used to find the indexes to test the attachment degree f at
    
    - 113-128: attachment degree f tests. The indexes to test at are found and then tested. A range of points was chosen, one pre alpha0, 2 in the linear region that should be fully attached (or starting to detach, depending), one at the cusp of deep stall, and one in deep stall. There is no relative tolerances for these tests, and there likely should be as webplot digitizer is a very inexact science. It would only need to be applied at places where it is not expected to be 1.0 or 0.0 (so 10deg and 18deg)-> currently both the 10deg and 18deg are failing, which is what the plot shows as well
    
    - 129-135: Critical alpha tests. I test, with a relative tolerance (currently 1%) for alphasep, alpha0 and dcldalpha and dcndalpha. It should be noted that these all pass, but that is because I used this file to get the values for them as I don't know where else to get the values, but it should work for future tests if something breaks. although the upper alphasep is hardcoded, as well as dcl and dcn, so those won't fail (See line 49).
    
    - 137-143: Setup tests. These are tests to make sure that I set everything up right and am using the correct separation function and solver etc. These aren't necessary and may want to be discarded once the other solvers are fully functional etc. For now it is a good way to minimize possible sources of error and discrepency.
  
  - [TestingTestTesting.jl](test\TestingTestTesting.jl)
    
    - This file was made to learn more about the Julia Test package. It goes through a bunch of different styles of tests and is heavily commented. It is not needed but could be helpful for me or someone in the future to study if help is needed creating or analyzing a test file.

- **testing/Oye** (ie the Oye folder in the testing folder)
  
  - [OyeComparer.jl](testing\Oye\OyeComparer.jl)
    
    - This file was created to compare the DSM model to Larsen's plot for the Vertol airfoil undergoing dynamic stall. It is concerning that it does not match up, and more digging needs to be done to see why not. The current output plot can be seen [here](Outputs\DSMvLarsenApr262023.png). There is a lag and the shape seems to be correct, but it is very far off of the actual values, more needs to be done checking the different methods. Messing with Tau (ie increasing the A input from .07 to like .2) makes the left hand side of the plot and part of the bottom match better, but the right hand side overshoots by a lot.
    
    - The code is well commented, but for more explanation-
    
    - 21-41: setting up  the airfoil and Oye structs
    
    - 42-60: Creating the angle of attack function (from Larsen's paper). This could be a source of the error and should be double checked
    
    - 62-63: Solving the model
    
    - 65-127: Cleaning up and plotting. I extract the polar from some csv files and put it in a format that is easy to plot
  
  - [SeparationFunctionNotes.md](testing\Oye\SeperationFunctionNotes.md)
    
    - This contains lots of notes that I took as I went through the airfoil file and tried to understand. Of most importance is the **Will Do List** at the top. It will be helpful in deciding what steps to take next! It has some good items and also provides a good history of what I did.
  
  - [OyeExample.jl](testing\Oye\OyeExample.jl)
    
    - This is outdated and was a precursor to OyeComparer.jl . I made it as I was following a tutorial
  
  - Outputs (folder)
    
    - This folder contains several critical outputs (and inputs as well unfortunately) for the various files, not just in the testing folder, but also elsewhere. Here are a few of note
    
    - [AllModelEquations.jpg](testing\Oye\Outputs\AllModelEquations.jpg) - This is a picture of the Larsen, Faber, and Hansen models and all their relevant equations written on a white board and is very helpful
    
    - [LarsenFig4cFdynAttachment.csv](testing\Oye\Outputs\LarsenFig4cFdynAttachment.csv) - This is the csv file from webplot digitizing the Larsen fdyn plot
    
    - [ClTypesPlotted.png](testing\Oye\Outputs\ClTypesPlotted.png) - The plots of all of the different Cl types (inviscid, fully separated etc) as well as the attachment degree from Larsen and DSM. It can be helpful in better understanding what is going on.
    
    - [LarsenFig4bStaticClPlotNaca63418.png](testing\Oye\Outputs\LarsenFig4bStaticClPlotNaca63418.png) - This is the Static Cl plot and will be important for future testing and comparing. It is already webplot digitized and saved in Polars as a csv. My test file uses it.
