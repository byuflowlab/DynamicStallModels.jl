# Weekly Update

### Mar 7th, 2023

- Fixed f (dynamic attachment degree) - an issue in the reverse separation point function 
  
  ![](C:\Users\child\.julia\dev\DynamicStallModels\testing\Oye\Outputs\fdynattachmentfix.png)

- Weston's fix helped me, but I am comparing inconsistent methods (Larsen's to Hansen's I think), which means the coefficients I am feeding in are wrong.

![image.png](C:\Users\child\.julia\dev\DynamicStallModels\testing\Oye\Outputs\updatedmodel.png)

- TestingTestTesting :) possibly the best file name I have ever come up with, I learned how to do several different types of tests and will be implementing it in a separation point test file
  
  ![](C:\Users\child\.julia\dev\DynamicStallModels\testing\Oye\Outputs\TestingTestTestingPic.png)

## March 15th, 2023

![](C:\Users\child\AppData\Roaming\marktext\images\2023-03-15-13-15-08-image.png)

Broke 2 sets of blades this morning! The first set I have seen!

![](C:\Users\child\.julia\dev\DynamicStallModels\testing\Oye\Outputs\DSMABitBetter.png)

We actually have a time lag and similar shape! I will start comparing f values and static curve values to debug-> good lesson on where to start

<img src="file:///C:/Users/child/.julia/dev/DynamicStallModels/testing/Oye/Outputs/FunctionCallMap.png" title="" alt="" width="655">

Function call tracker, so I could make sure everything was going where I thought it was and that I replaced every hardcoded separation point equation with a function call

Created To-Do List/Plan and made progress! 

### March 29, 2023

![](C:\Users\child\AppData\Roaming\marktext\images\2023-03-29-12-45-19-image.png)

Blue is Larsen's data, orange is our data we are overshooting for some reason

![](C:\Users\child\AppData\Roaming\marktext\images\2023-03-29-12-46-57-image.png)

I wrote a test file, there are still some bugs I need to work out, it caught the issue in the above plot, but not as well as I would have liked, so I may need to add some more test cases. I'm not quite sure why the critical alphas are failing and will need to look into that, It is likely I just set it up wrong or am comparing it to an incorrect value

- Accomplished a big chunk of my To Do list
- Learned to document everything! (I didn't document where my critical alpha test values came from, and I had to spend a lot of time relearning how to plot f values and figuring out why it didn't match my old plots, I was comparing different methods)

## April 12th, 2023

![20230408_002553.jpg](C:\Users\child\.julia\dev\DynamicStallModels\testing\Oye\Outputs\ClTypesHandDraw.jpg)

![20230408_002553.jpg](C:\Users\child\.julia\dev\DynamicStallModels\testing\Oye\Outputs\ClTypesPlotted.png)
