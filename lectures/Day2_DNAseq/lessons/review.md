---

title: "Review of Day 1"
author: "Shahin Shahsavari"
date: March 2020

---

# Day 1 Review

Yesterday, we learned the following commands to get around in the linux terminal:

```bash
pwd
cd
ls
ls -l -h
```

We then introduced the following for viewing contents of our files:

```bash
cat
head
tail
head -n 20
```

And lastly, if you want to move, copy, rename, create, or remove files and directories, we learned:

```bash
cp
cp -r # to copy folders
mdkir
mv
rm
rm -r # to remove folders
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
# Send the output of the previous command into a new program aka pipe
`|`

```

We also had

```bash
wc -l
cut -f 1,3 # extract first and third columns
```
