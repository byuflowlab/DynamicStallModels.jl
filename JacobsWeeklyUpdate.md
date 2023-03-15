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
