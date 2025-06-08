
# Contributing to microbetag

Thank you for your interest in contributing to microbetag!
This guide outlines the recommended workflow for contributing code or documentation.

## 🚀 Getting Started

1. **Fork the Repository**

Click the **Fork** button on the top right of [microbetag](https://github.com/msysbio/microbetag) 
and **clone your fork locally**:

    git clone https://github.com/<your-username>/microbetag.git
    cd microbetag

2. Set Up Remotes

Add the upstream repository to keep your fork in sync:

    git remote add upstream https://github.com/msysbio/microbetag.git
    git remote -v

This should show:

    origin    https://github.com/<your-username>/microbetag.git (fetch)
    origin    https://github.com/<your-username>/microbetag.git (push)
    upstream  https://github.com/msysbio/microbetag.git (fetch)
    upstream  https://github.com/msysbio/microbetag.git (push)

## 🔄 Syncing with the Upstream Repository

Before starting a new feature or fix, make sure your local repository is up to date:

    git fetch upstream
    git checkout develop
    git pull upstream develop --rebase

## 🌱 Creating a Feature Branch

Use a descriptive name for your branch:

    git checkout -b feature/my_descriptive_branch_name

Alternatively, if you're continuing from a prior state:

    git branch feature/the_fastest_sampling_algo_ever
    git checkout feature/the_fastest_sampling_algo_ever

## 🚀 Pushing to Your Fork

Push your branch to your fork (replace origin with your fork name if different):

    git push -u origin feature/my_descriptive_branch_name

The `-u` flag sets the upstream tracking branch, so future git push/pull commands are simpler.


## ✅ Submitting a Pull Request (PR)

   1. Go to your fork on GitHub.

   2. Click **Compare & pull request**.

   3. Base branch should be `develop`, and compare should be your feature branch.

   4. Add a clear title and description explaining your changes.

   5. Link any related issues or discussions.

   6. Submit the pull request.

## 🔧 Tips

 1. Write clear commit messages.

 2. Test your code before opening a PR; it's a good thing to provide a new test when building a new feature.

 3. Use draft PRs if you're still working but want early feedback.

