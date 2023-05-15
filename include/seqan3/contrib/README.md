The `std` directory is a git subtree and includes the content of https://github.com/seqan/seqan-std.

For a general overview regarding `git subtree`, see, for example, https://www.atlassian.com/git/tutorials/git-subtree.

Pulling upstream changes can be achieved by
```
git subtree pull --prefix include/seqan3/contrib/std https://github.com/seqan/seqan-std.git main --squash
```

Similarly, you can push subtree changes to your fork:
```
git subtree push --prefix include/seqan3/contrib/std https://github.com/<user>/seqan-std.git <branch>
```
