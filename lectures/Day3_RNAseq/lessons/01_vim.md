---
title: Introduction to Vim
author: "Shahin Shahsavari"
date: "November 2021"
---

## Learning Objectives

* Learn basic operations using the Vim text editor

## Writing files

We've been able to do a lot of work with files that already exist, but what if we want to write our own files. Obviously, we're not going to type in a FASTA file, but you'll see as we go, there are a lot of reasons we'll want to write/create a file or edit an existing file.

To create or edit files we will need to use a **text editor**. When we say, "text editor," we really do mean "text": these editors can
only work with plain character data, not tables, images, or any other
media. The types of text editors available can generally be grouped into **graphical user interface (GUI) text editors** and **command-line editors**.

### GUI text editors

A GUI is an interface that has buttons and menus that you can click on to issue commands to the computer and you can move about the interface just by pointing and clicking. You might have heard the names of GUI (graphical) text editors, such as [TextWrangler](http://www.barebones.com/products/textwrangler/), [VSCode](https://code.visualstudio.com/), and [Notepad++](http://notepad-plus-plus.org/), which allow you to write and edit plain text documents. These editors often have features to easily search text, extract text, and highlight syntax from multiple programming languages. They are great tools, but since they are 'point-and-click', we cannot efficiently use them from the command line remotely on a compute cluster.

### Command-line editors

When working remotely, we need a text editor that functions from the command line interface. Within these editors, since you cannot 'point-and-click', you must navigate the interface using the arrow keys and shortcuts. 

While there are simpler editors available for use (i.e. [nano](http://www.nano-editor.org/)), most computational scientists tend to favor editors that have greater functionality. Some popular editors include [Emacs](http://www.gnu.org/software/emacs/), [Vim](http://www.vim.org/), or a graphical editor such as [Gedit](http://projects.gnome.org/gedit/). These are editors which are generally available for use on high-performance compute clusters.

### Introduction to Vim 

To write and edit files, we're going to use a text editor called 'Vim'. Vim is a very powerful text editor, and it offers extensive text editing options. However, in this introduction we are going to focus on exploring some of the more basic functions. There is a lot of functionality that we are not going to cover during this session, but encourage you to go further as you become more comfortable using it. To help you remember some of the keyboard shortcuts that are introduced below and to allow you to explore additional functionality on your own, there is a great [cheatsheet](https://vim.rtorr.com/). Please rest assured that **most** bioinformaticians and computer scientists find it hard to use vim, though it is worth learning if you spend lots of time developing on remote servers.


### Vim Interface

You can create a document by calling a text editor and providing the name of the document you wish to create. Change directories to the `unix_lesson/other` folder and create a document using `vim` entitled `draft.txt`:

```bash
$ cd ~/Day3/text_editor
	
$ vim draft.txt
```

Notice the `"draft.txt" [New File]` typed at the bottom left-hand section of the screen. This tells you that you just created a new file in vim. 


### Vim Modes
Vim has **_two basic modes_** that will allow you to create documents and edit your text:   

- **_command mode (default mode):_** will allow you to save and quit the program (and execute other more advanced commands).  

- **_insert (or edit) mode:_** will allow you to write and edit text


Upon creation of a file, vim is automatically in command mode. Let's _change to insert mode_ by typing <kbd>i</kbd>. Notice the `--INSERT--` at the bottom left hand of the screen. Now type in a few lines of text:

<img src="../img/vim_insert.png" width="600">

After you have finished typing, press <kbd>esc</kbd> to enter command mode. Notice the `--INSERT--` disappeared from the bottom of the screen.

### Vim Saving and Quitting
To **write to file (save)**, type <kbd>:w</kbd>. You can see the commands you type in the bottom left-hand corner of the screen. 

<img src="../img/vim_save.png" width="600">

After you have saved the file, the total number of lines and characters in the file will print out at the bottom left-hand section of the screen.

<img src="../img/vim_postsave.png" width="600">

Alternatively, we can **write to file (save) and quit**. Let's do that by typing <kbd>:wq</kbd>. Now, you should have exited vim and returned back to your terminal window.

To edit your `draft.txt` document, open up the file again by calling vim and entering the file name: `vim draft.txt`. Change to insert mode and type a few more lines (you can move around the lines using the arrows on the keyboard). This time we decide to **quit without saving** by typing <kbd>:q!</kbd>
 
<img src="../img/vim_quit.png" width="600">


While we cannot point and click to navigate the document, we can use the arrow keys to move around. Navigating with arrow keys can be very slow, so Vim has shortcuts (which are completely unituitive, but very useful as you get used to them over time). Check to see what mode you are currently in. While in command mode, try moving around the screen and familarizing yourself with some of these shortcuts:

| key              | action                 |
| ---------------- | ---------------------- |
| <button>gg</button>     | to move to top of file |
| <button>G</button>     | to move to bottom of file     |
| <button>$</button>     | to move to end of line |
| <button>0</button>     | to move to beginning of line     |
| <button>w</button>     | to move to next word     |
| <button>b</button>     | to move to previous word     |


In addition to shortcuts for navigation, vim also offers editing shortcuts such as:

| key              | action                 |
| ---------------- | ---------------------- |
| <button>dw</button>     | to delete word |
| <button>dd</button>     | to delete line     |
| <button>u</button>     | to undo |
| <button>Ctrl + r</button>     | to redo     |
| <button>/*pattern*</button>     | to search for a pattern (*n/N* to move to next/previous match)    |
| <button>:%s/*search*/*replace*/g</button>     | to search for a pattern and replace for all occurences     |

Practice some of the editing shortcuts, then quit the document without saving any changes.

*** 

**Exercise**

We have covered some basic commands in vim, but practice is key for getting comfortable with the program. Let's
practice what we just learned in a brief challenge.

1. Navigate into your "~/Day1/R_data"
2. Use `cat` to view the content of patients.tsv
3. Use vim to open the file and then close without saving
4. Open the file again and change *harvey*'s height from NA to 174 
5. Save the file.
6. Open the file again and add one more person to the end of your file
	- We will read this file in R later

***

### Overview of vim commands

**Vim modes:**

| key              | action                 |
| ---------------- | ---------------------- |
| <button>i</button>     | insert mode - to write and edit text |
| <button>esc</button>     | command mode - to issue commands / shortcuts  |


**Saving and quiting:**

| key              | action                 |
| ---------------- | ---------------------- |
| <button>:w</button>     | to write to file (save) |
| <button>:wq</button>     | to write to file and quit     |
| <button>:q!</button>     | to quit without saving |


**Shortcuts for navigation:**

| key              | action                 |
| ---------------- | ---------------------- |
| <button>gg</button>     | to move to top of file |
| <button>G</button>     | to move to bottom of file     |
| <button>$</button>     | to move to end of line |
| <button>0</button>     | to move to beginning of line     |
| <button>w</button>     | to move to next word     |
| <button>b</button>     | to move to previous word     |

**Shortcuts for editing:**

| key              | action                 |
| ---------------- | ---------------------- |
| <button>dw</button>     | to delete word |
| <button>dd</button>     | to delete line     |
| <button>u</button>     | to undo |
| <button>Ctrl + r</button>     | to redo     |
| <button>:set number</button>     | to number lines |
| <button>:set nonumber</button>     | to remove line numbers    |
| <button>/pattern</button>     | to search for a pattern (*n/N* to move to next/previous match)    |
| <button>:%s/search/replace/g</button>     | to search for a pattern and replace for all occurences     |	

***
