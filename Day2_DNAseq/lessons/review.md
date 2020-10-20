Title: "Review of Day 1"
Date: October 2020

# Day 1 Review

Yesterday, we learned the following commands to get around in the linux terminal:

```bash
pwd
cd
ls
ls -lh
ls -la
```

We then introduced the following for viewing contents of our files:

```bash
cat
less
head
tail
```

And lastly, if you want to move, copy, rename, create, or remove files and directories, we learned:

```bash
cp
cp -r
mdkir
mv
rm
```

For searching through files, and redirection:

```bash
grep
# two line after a string match and one line before
grep -A 2 -B 1

# Redirect the terminal output and create or rewrite a file
`>`
# Append
`>>`
# Send the output of a command into a new program
`|`

```

We also had

```bash
wc -l
cut -f 1,3
sort
```


We also learned vim


**Practice**

1. In your ~/unix\_lesson folder, use vim to create a shell script "practice.sh"
2. we want this shell script to do the following:
- store "/data/review" into a variable called "dir"
- use the dir variable to list all the files in that directory that end with ".tsv"
- do the same for files that end with ".csv"


