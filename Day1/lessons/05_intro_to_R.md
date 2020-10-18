---
title: "Intro to R"
author: "Shahin Shahsavari"
date: "October 2020"
---

## Learning Objectives
- 
- How do you use it? 
  - Getting around the Unix file system
  - looking at files
  - manipulating files
  - automating tasks
- What is it good for?

## 1. What is Linux

Linux is a family of open source Unix-like operating systems based on the Linux kernel, an operating system kernel first released in September, 1991, by Linus Torvalds. 67% of the world’s servers run some variant of Unix or Linux. Android also uses a modified version of the Linux Kernel. There are a wide variety of Linux distributions (Ubuntu, Mint, Fedora, Arch, etc).

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
- 2020: All 500 of the top 500 supercomputers in the world run variants of Linux today.


## 2. Setting up

To further explore Linux, we are going to log into our teaching server that runs Ubuntu Linux.

### Logging in using X2Go Client

1. open the **X2Go Client**
2. click on session -> New Session
3. Enter the following:
        * Host: 3.237.4.148
        * Login: USERNAME
        * Session type: XFCE
