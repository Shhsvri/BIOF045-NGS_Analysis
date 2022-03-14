---
title: "The Shell"
author: "Shahin Shahsavari"
date: "March 2022"
---

## Learning Objectives
- How do you access the shell/terminal?
- How do you use it?
  - getting around the Unix file system
  - viewing the content of files
  - manipulating files
  - automating tasks
- What is it good for?


## 1. What is Linux

Linux is a family of open source Unix-like operating systems based on the Linux kernel, an operating system kernel
first released in September, 1991, by Linus Torvalds. 67% of the world’s servers run some variant of Unix or Linux.
Android also uses a modified version of the Linux Kernel. There are a wide variety of Linux distributions (Ubuntu,
Mint, Fedora, Arch, etc).


### 1.1 A Brief History of Unix

- 1969: Unix was developed at Bell Labs by Ken Thompson, Dennis Ritchie & others.
- 1971: First Edition used for text processing of patent documents.
- 1973: Ritchie and Thompson developed C, with which they wrote the Unix operating system.
- 1983: System V is the first supported release by AT&T with an installation base of 45,000.
- 1991: Linus Torvalds commences Linux Development. Solaris 1.0 debuts.
- 1993: Red Hat Linux and Debian Linux were released which are widely used.
- 1997: Titanic is the first major film to be largely produced on Linux servers.
- 1999: Linux Kernel 2.2 released to celebrate Unix's 30th birthday.
- 2001: Apple releases Mac OS X, based on BSD Unix.
- 2003: Linux 2.6 kernel released. Red Hat creates RHEL and Fedora.
- 2004: Ubuntu, a popular desktop Linux distribution based on Debian, first released.
- 2007: Mac OS X certified to Unix standard.
- 2010: Apple reports 50 million desktops and growing, all Certified Unix Systems.
- 2012: 500 million Android (Linux) and iOS (Darwin) devices have been sold world wide.
- 2020: All 500 of the top 500 supercomputers in the world run variants of Linux.


## 2. Setting up

To further explore Linux, we are going to log into our teaching server that runs Ubuntu Linux.


### Logging in using X2Go Client

1. open the **X2Go Client**
2. click on session -> New Session
3. Enter the following:
	* Host: 54.234.225.44
	* Login: USERNAME
	* Session type: XFCE

<img src="../img/X2Go_login.png" width="600">

Click on the session you created, enter password, and log in.

You should see the content of your Desktop in a new window.


## 3. Getting started with the Command Line (Terminal/Bash Shell)

In the *Applications* drop down menu, open the terminal emulator.

<img src="../img/terminal_emulator.png" width="600">

This new window is referred to as the terminal or the shell. If you learn some basic commands, you can get some complicated jobs done with them.

The command prompt will have some characters before it, something like `shahin@ip-172-31-75-209:~$ `, this is telling what your username and the internal address of the cloud instance you are working on is.

Let's try a few commands:

```bash
$ date
$ cal
$ whoami
```


## 3.1 Unix Filesystem

Just like a file manager, you can navigate through your computers and move around in folders (directories) using
the terminal. Let's see where we are:

```bash
$ pwd
```

> 'pwd' stands for 'print working directory'. It tells you what folder you are currently located in.

Let's list all the files and folders we have here.

```bash
$ ls
```

> ls stands for 'list' and it lists the contents of a directory.

Let's go into our `Day1` directory and see what we copied there. Type:

```bash
$ cd Day1

$ ls
```

You will see:

```
README.txt  R_data  genomics_data  hello  other  raw_fastq  reference_data
```

There are seven items listed. What types of files are they? We can use a "modifier" with `ls` to get more information; this modifier is called an argument (more below).

```bash
$ ls -F

README.txt  R_data/  genomics_data/  hello*  other/  raw_fastq/  reference_data/
```

Anything with a "/" after it is a directory. If there are no decorations after the name, it's a file. Files that end with "\*" are executables.

> All commands are essentially programs that are able to perform specific, commonly-used tasks.

You can also use the command:

```bash
$ ls -l
```

to see whether items in a directory are files or directories. `ls -l` gives a lot more information too.
```
total 44
-rw-r--r-- 1 shahin shahin 379 Jul 11 11:31 README.txt
drwxr-xr-x 2 shahin shahin 4096 Jul 11 11:31 R_data
drwxr-xr-x 2 shahin shahin 4096 Jul 11 11:31 genomics_data
-rwxr-xr-x 1 shahin shahin 16840 Jul 11 11:31 hello
drwxr-xr-x 2 shahin shahin 4096 Jul 11 11:31 other
drwxr-xr-x 2 shahin shahin 4096 Jul 11 11:31 raw_fastq
drwxr-xr-x 2 shahin shahin 4096 Jul 11 11:31 reference_data
```

Let's go into the raw\_fastq directory and see what is in there.

```bash
$ cd raw_fastq/

$ ls -F

treated_1.fastq  treated_3.fastq    untreated_2.fastq
treated_2.fastq  untreated_1.fastq  untreated_3.fastq
```

The six items in this directory have no trailing slashes, so they are all files, not folders or programs.


#### Arguments

Most commands take additional arguments that control their exact behavior. For example, `-F` and `-l` are arguments to `ls`. The `ls` command, like many commands, take a lot of arguments. Another useful one is `-a`, which shows everything, including hidden files.  How do we know what the available arguments that go with a particular command are?

Most commonly used shell commands have a manual available in the shell. You can access the
manual using the `man` command. Try entering:

```bash
$ man ls
```

This will open the manual page for `ls`. Use the 'space' key to go forward and 'b' to go backwards. When you are done reading, just hit `q` to quit.

Commands that are run from the shell can get extremely complicated. To see an example, open up the manual page for the `find` command. No one can possibly learn all of these arguments, of course. So you will probably find yourself referring to the manual page frequently.

> If the manual page within the terminal is hard to read and traverse, the manual exists online,
	 use your web searching powers to get it! In addition to the arguments, you can also find good examples online; Google is your friend.


## The Unix directory file structure
 
As you've already just seen, you can move around in different directories or folders at the command line. Why would you want to do this, rather than just navigating around the normal way using a GUI (GUI = Graphical User Interface, pronounced like "gooey")?

#### Moving around the file system

Let's practice moving around a bit.

We're going to work in that `Day1` directory.

First we did something like go to the folder of our username. Then we opened `Day1` then `raw_fastq`

Like on any computer you have used before the file structure within unix is hierarchical, like an upside down tree with root (/) as the starting point of the tree-like structure:

<img src="../img/directory_structure.png" width="600">

That root (/) is often also called the 'top' level.

When you log in to a remote computer you are on one of the branches of that tree, your home directory (e.g. /home/USERNAME)

> On mac OS, which is a UNIX-based OS, the root level is also "/". On a windows OS, it is drive specific; generally "C:\" is considered root, but it changes to "D:/", if you are on that drive.

Now let's go do that same navigation at the command line.

Type:

```bash
$ cd ~
```

> No matter where you are in the directory system, `cd ~` will always bring you back to your home directory. `cd /home/USERNAME`

---

**Exercise 1**

Now using `cd` and `ls`, go in to the `Day1` directory and list its contents. From there, go into the `raw_fastq` directory, and list its contents.

---

Let's also check to see where we are. Sometimes when we're wandering around in the file system, it's easy to lose track of the directory structure. The command that tells you this is:

```bash
$ pwd
```


What if we want to move back up and out of the `raw_fastq` directory? Can we just type `cd Day1`? Try it and see what happens.

To go 'back up a level' we can use `..`

Type:

```bash
$ cd ..
```

Use `pwd` to fidn out where you are and `ls` to list the content of the current directory. 

> `..` denotes parent directory, and you can use it anywhere in the system to go back to the parent directory. Can you think of an example when this won't work?

Finally, there is handy command that can help you see the structure of any directory, namely `tree`.

```bash
# Ensure that you are in your Day1 directory, then run the following command

$ tree
```

#### Examining the contents of other directories

By default, the `ls` commands lists the contents of the working directory (i.e. the directory you are in). You can always find the directory you are in using the `pwd` command. However, you can also give `ls` the names of other directories to view. Navigate to the home directory if you are not already there.

Type:

Then enter the command:

```bash
$ ls Day1/
```

This will list the contents of the `Day1` directory without you having to navigate there.

The `cd` command works in a similar way.

```bash
$ cd Day1/raw_fastq/
$ pwd
```

You should now be in `raw_fastq` and you got there without having to go through the intermediate directory. 

> If you are aware of the directory structure, you can string together as long a list as you like.


---

**Exercise 2**

List the `treated_1.fastq` file from your home directory without changing directories

---

## Full vs. Relative Paths

The `cd` command takes an argument which is the directory name. Directories can be specified using either a *relative path* or a *full path*. As we know, the directories on the computer are arranged into a hierarchy. The full path tells you where a directory is in that hierarchy. Navigate to the home directory (`cd`). Now, enter the `pwd` command and you should see:

```bash
$ pwd
```

```
/home/USERNAME
```

	which is the full path for your home directory. This tells you that you are in a directory called `username`, which sits inside a directory called `home` which sits inside the very top directory in the hierarchy, the *root directory*. So, to summarize: `username` is a directory in `home` which is a directory in `/`.

	Now enter the following command:

	```bash
	$ cd /home/USERNAME/Day1/raw_fastq/
	```

	This jumps to `raw_fastq`. Now go back to the home directory (`cd ~`). We saw
	earlier that the command:

	```bash
	$ cd Day1/raw_fastq/
	```

	had the same effect - it took us to the `raw_fastq` directory. But, instead of specifying the full path (`/home/USERNAME/Day1/Day1/raw_fastq`), we specified a *relative path*. In other words, we specified the path **relative to our current working directory**. 

	**A full path always starts with a `/`, a relative path does not.**

	A relative path is like getting directions from someone on the street. They tell you to "go right at the Stop sign, and then turn left on Main Street". That works great if you're standing there together, but not so well if you're trying to tell someone how to get there from another country. A full path is like GPS coordinates. It tells you exactly where something is no matter where you are right now.

	You can usually use either a full path or a relative path depending on what is most convenient. If we are in the home directory, it is more convenient to just enter the relative path since it involves less typing.

	Over time, it will become easier for you to keep a mental note of the structure of the directories that you are using and how to quickly navigate among them.

	---

	**Exercise 3**

	Change directories to `/home/USERNAME/Day1/raw_fastq/`, and list the contents of `Day1/other` without changing directories again.

---

### Saving time with tab completion, wildcards and other shortcuts 

#### Tab completion

Navigate to the home directory. Typing out directory names can waste a lot of time. When you start typing out the name of a directory, then hit the tab key, the shell will try to fill in the rest of the directory name. For example, type `cd` to get back to your home directly, then enter:

```bash
$ cd Da<tab>
```

The shell will fill in the rest of the directory name for `Day1`. Now go to `Day1/raw_fastq` and 

```bash
$ ls treated<tab><tab>
```

When you hit the first tab, nothing happens. The reason is that there are multiple directories in the home directory which start with `treated_`. Thus, the shell does not know which one to fill in. When you hit tab again, the shell will list the possible choices.

Tab completion can also fill in the names of commands. For example, enter `e<tab><tab>`. You will see the name of every command that starts with an `e`. One of those is `echo`. If you enter `ech<tab>` you will see that tab completion works. 

> **Tab completion is your friend!** It helps prevent spelling mistakes, and speeds up the process of typing in the full command.

#### Wild cards

Navigate to the `~/Day1/raw_fastq` directory. This directory contains FASTQ files from a next-generation sequencing dataset. 

The '*' character is a shortcut for "everything". Thus, if you enter `ls *`, you will see all of the contents of a given directory. Now try this command:

```bash
$ ls treated*
```

This lists every file that ends with a `fastq`. This command:

```bash
$ ls /usr/bin/*.sh
```

Lists every file in `/usr/bin` that ends in the characters `.sh`.

```bash
$ ls treated*.fastq
```

> lists only the files that begin with 'treated' and end with 'fastq'

So how does this actually work? The shell (bash) considers an asterisk "*" to be a wildcard character that can be used to substitute for any other single character or a string of characters. 

> An asterisk/star is only one of the many wildcards in UNIX, but this is the most powerful one and we will be using this one the most for our exercises.

****

**Exercise**

Do each of the following using a single `ls` command without
navigating to a different directory.

1.  List all of the files in `/bin` that start with the letter 'c'
2.  List all of the files in `/bin` that contain the letter 'a'
3.  List all of the files in `/bin` that end with the letter 'o'

BONUS: List all of the files in `/bin` that contain the letter 'a' or 'c'.

****


#### Shortcuts

There are some shortcuts which you should know about. Dealing with the
home directory is very common. So, in the shell the tilde character,
"~", is a shortcut for your home directory. Navigate to the `raw_fastq`
directory:

```bash
$ cd
```

```bash
$ cd Day1/raw_fastq
```

Then enter the command:

```bash
$ ls ~
```

This prints the contents of your home directory, without you having to type the full path because the tilde "~" is equivalent to "/home/username".

Another shortcut is the "..":

```bash
$ ls ..
```

The shortcut `..` always refers to the directory above your current directory. So, it prints the contents of the `Day1`. You can chain these together, so:

```bash
$ ls ../..
```

prints the contents of `/home/username` which is your home directory. 

Finally, the special directory `.` always refers to your current directory. So, `ls`, `ls .`, and `ls ././././.` all do the same thing, they print the contents of the current directory. This may seem like a useless shortcut right now, but we used it earlier when we copied over the data to our home directory.


To summarize, while you are in your home directory, the commands `ls ~`, `ls ~/.`, and `ls /home/username` all do exactly the same thing. These shortcuts are not necessary, but they are really convenient!

#### Command History

You can easily access previous commands.  Hit the up arrow. Hit it again.  You can step backwards through your command history. The down arrow takes your forwards in the command history.

You can also review your recent commands with the `history` command.  Just enter:

```bash
$ history
```

to see a numbered list of recent commands, including this just issues
`history` command. Only a certain number of commands are stored and displayed with `history`, there is a way to modify this to store a different number.

> **NOTE:** So far we have only run very short commands that have few or no arguments, and so it would be faster to just retype it than to check the history. However, as you start to run analyses on the commadn-line you will find your commands to be more complex and the history to be very useful!

**Other handy command-related shortcuts**

- <button>Ctrl + C</button> will cancel the command you are writing, and give you a fresh prompt.
- <button>Ctrl + A</button> will bring you to the start of the command you are writing.
- <button>Ctrl + E</button> will bring you to the end of the command.

## Examining Files

We now know how to move around the file system and look at the
contents of directories, but how do we look at the contents of files?

The easiest way to examine a file is to just print out all of the
contents using the command `cat`. Print the contents of `Day1/other/sequences.fa` by entering the following command:

```bash
$ cat ~/Day1/other/sequences.fa
```

This prints out the all the contents of `sequences.fa` to the screen.

> `cat` stands for catenate; it has many uses and printing the contents of a files onto the terminal is one of them.

What does this file contain?

`cat` is a terrific command, but when the file is really big, it can be annoying to use.
In practice, when you are running your analyses on the command-line you will most likely be dealing with large files.
The command, `less`, is useful for this case. Let's take a look at the list of raw_fastq files and add the `-h` modifier:

```bash
ls -l -h ~/Day1/raw_fastq
```

> The `ls` command has a modifier `-h` when paired with `-l`, will print sizes of files in human readable format 

In the fourth column you see the size of each of these files, and you can see they are quite large, so we probably do not want to use the `cat` command to look at them. Instead, we can use the `less` command. 

Move back to our `raw_fastq` directory and enter the following command:

```bash
less treated_1.fastq
```

We will explore FASTQ files in more detail later, but notice that FASTQ files have four lines of data associated with every sequence read. Not only is there a header line and the nucleotide sequence, similar to a FASTA file, but FASTQ files also contain quality information for each nucleotide in the sequence. 


```bash
$ head treated_1.fastq
```

```bash
$ tail treated_1.fastq
```

The `-n` option to either of these commands can be used to print the first or last `n` lines of a file. To print the first/last line of the file use:

```bash
$ head -n 1 treated_1.fastq

$ tail -n 1 treated_1.fastq
```

## Creating, moving, copying, and removing

Now we can move around in the file structure, look at files, search files, redirect. But what if we want to do normal things like copy files or move them around or get rid of them. Sure we could do most of these things without the command line, but what fun would that be?! Besides it's often faster to do it at the command line, or you'll be on a remote server like Amazon where you won't have another option.

Our raw data in this case is fastq files. We don't want to change the original files, so let's make a copy to work with.

Lets copy the file using the copy `cp` command. Navigate to the `raw\_fastq` directory and enter:

```bash
$ cp treated_1.fastq treated_1-copy.fastq

$ ls -l
```

Now `treated_1-copy.fastq` has been created as a copy of `treated_1.fastq`

Let's make a 'backup' directory where we can put this file.

The `mkdir` command is used to make a directory. Just enter `mkdir`
followed by a space, then the directory name.

```bash
$ mkdir backup
```

> File/directory/program names with spaces in them do not work in unix, use characters like hyphens or underscores instead.

We can now move our backed up file in to this directory. We can move files around using the command `mv`. Enter this command:

```bash
$ mv *copy.fastq backup/
```

```bash
$ ls -l backup/

-rw-rw-r-- 1 shahin shahin 75706556 Oct 11 13:56 treated_1-copy.fastq
```

The `mv` command is also how you rename files. Since this file is so
important, let's rename it:

```bash
$ cd backup

$ mv treated_1-copy.fastq treated_1-backup.fastq

$ ls

treated_1-backup.fastq
```

Finally, we decided this was silly and want to start over.

```bash
$ cd ..

$ rm backup/treated*.fastq
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

## Commands, options, and keystrokes covered

```
bash
cd
ls
man
pwd
~           # home dir
.           # current dir
..          # parent dir
*           # wildcard
echo
ctrl + c    # cancel current command
ctrl + a    # start of line
ctrl + e    # end of line
history
cat
less
head
tail
cp
mkdir
mv
rm
```

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
