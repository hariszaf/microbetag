## Packaging

We use `poetry`. 

### For the development part

Run `poetry shell` to activate the env. 

Once a depedency is added use `poetry add` to add the new package in your `.toml`.

To build the `microbetag` package you may run:

```
poetry build 
```

and then, if you want to make it callable from anywhere on your pc:

```
pip install dist/*.tar.gz 
```
