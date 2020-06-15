git branch -d binder --force
git checkout -b binder
git branch -u origin/binder
echo "# Python dependencies for mybinder
dependencies:
- matplotlib
- numpy
- pip
- pip:
  - julia" >> environment.yml
julia --project binder_dependencies.jl
git add -A
git commit -m "added binder dependencies" -a
git push --force
rm environment.yml
git checkout master
