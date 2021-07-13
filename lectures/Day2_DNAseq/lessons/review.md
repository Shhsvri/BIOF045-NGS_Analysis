---

title: "Review of Day 1"
author: "Shahin Shahsavari"
date: July 2020

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

---

**Exercise**

1.  Change your directory to `~/Day2/Day1_review`. `cd`
2.  Find the size of all the files in the `Day1_review` folder. `ls`
3.  Find the number of lines in `ptA_R1.fastq`. `wc`
4.  How many reads in `ptA_R1.fastq` contain `NNNNN`? `grep`, `|` and `wc`
5.  Write the fastq records that contain `NNNNN` into a new file called `bad_reads.fastq`.

---

If you want to move, copy, rename, create, or remove files and directories, we use:

```bash
cp
cp -r # to copy folders
mdkir
mv
rm
rm -r # to remove folders
```

Let's say we decide to create a folder called `backup` and move a copy of our fastq files there.

First let's create a copy of our original fastq file:

```
cp ptA_R1.fastq ptA_R1-copy.fastq
cp bad_reads.fastq bad_reads-copy.fastq
```

Then we create a backup folder. The `mkdir` command is used to make a directory. Just enter `mkdir` followed by a space, then the directory name that you would like to create.
```
$ mkdir backup
```

> File/directory/program names with spaces in them do not work in unix, use characters like hyphens or underscores instead.

We can now move our backed up file in to this directory. We can move files around using the command `mv`. Enter this command:

```
$ mv *copy.fastq backup
$ ls -l backup
```

Finally, we decided that we are running out of space on our server and need to remove some files.

```
rm backup/*copy.fastq
```

> The `rm` file permanently removes the file. Be careful with this command. The shell doesn't
just nicely put the files in the Trash. They're really gone.
>
> Same with moving and renaming files. It will **not** ask you if you are sure that you want to "replace existing file". You can use `rm -i` if you want it to ask before deleting the file(s).

We really don't need these backup directories, so, let's delete both. By default, `rm`, will NOT delete directories, but you use the `-r` flag if you are sure that you want to delete the directories and everything within them.

```bash
$ rm -r backup/
```

- `-r`: recursive, commonly used as an option when working with directories, e.g. with `cp`.


#### Information on the shell

shell cheat sheets:<br>
* [http://fosswire.com/post/2007/08/unixlinux-command-cheat-sheet/](http://fosswire.com/post/2007/08/unixlinux-command-cheat-sheet/)
* [https://github.com/swcarpentry/boot-camps/blob/master/shell/shell_cheatsheet.md](https://github.com/swcarpentry/boot-camps/blob/master/shell/shell_cheatsheet.md)

Explain shell - a web site where you can see what the different components of
a shell command are doing.
* [http://explainshell.com](http://explainshell.com)
* [http://www.commandlinefu.com](http://www.commandlinefu.com)

Data Carpentry tutorial: [Introduction to the Command Line for Genomics](https://datacarpentry.org/shell-genomics/)

General help:
- http://tldp.org/HOWTO/Bash-Prog-Intro-HOWTO.html
- man bash
- Google - if you don't know how to do something, try Googling it. Other people
have probably had the same question.
- Learn by doing. There's no real other way to learn this than by trying it
out.  Write your next paper in vim (really emacs or vi), open pdfs from
the command line, automate something you don't really need to automate.

---

*This lesson has been developed by Shahin Shahsavari, adapted from [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

* *The materials used in this lesson were derived from work that is Copyright © Data Carpentry (http://datacarpentry.org/).
All Data Carpentry instructional material is made available under the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0).*
