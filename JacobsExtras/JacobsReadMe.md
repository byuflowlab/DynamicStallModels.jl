JacobsReadMe.md
**Important Local Locations:**
*Assuming I am currently working in FlowLab_DynamicStall\DynamicStallModels.jl*

1. DynamicStallModels package (this is the package used/referenced every time I run)
    [here](C:\Users\child\.julia\dev\DynamicStallModels)

2. The things I am currently working on will all be located in the same place I assume I am currently working.
    [here](C:\Users\child\Documents\Projects\FlowLab_DynamicStall\DynamicStallModels.jl)

**Plan:**
I will run verification tests and learn while in [2.](C:\Users\child\Documents\Projects\FlowLab_DynamicStall\DynamicStallModels.jl), and then once I am ready to make actual changes I will change to [1.](C:\Users\child\.julia\dev\DynamicStallModels) and change the actual package, pushing and pulling from there. 
*Note: The changes that matter most (in the actual package) are under the [src]("C:\Users\child\.julia\dev\DynamicStallModels\src") folder, but as mentioned above, that should not happen until I am ready.*

*Also note: the FlowLab_DynamicStall folder **is not** a github repo or branch etc, it is **only saved locally**. Only the DynamicStallModels.jl folder is, that means I must cd to that folder for the purposes of git, and make sure that I don't save anything locally that I don't need to.

**Github Notes:**
I have permission, but *must not* push to master and merge a pull request with master, unless explicitly cleared with Adam. If my branch is behind Master and I want to update mine I take the following steps:
1. `git checkout master` (this switches me to the master branch)
2. `git pull` (this updates my local master branch, this is necessary as I am operating locally and must update master before updating mine)
3. `git checkout jacobchild` (this switches me back to my branch)
4. `git merge master` (this merges master into my branch)
5. `git push` (pushes my merge from local to online) 

*further notes:*
1. To render markdown in VSCode, do ctrl+k followed quickly by v
2. @JacobChild does this give me some sort of notification?