#!/bin/sh

# changes to wiki/ directory
if [ -d ../../wiki ]; then
  cd ../../wiki/
else
  echo "wiki directory does not exist, nothing to do..."
  exit 1
fi

# checks if translation file exists
if [ ! -f translate_user_manual_to_markdown.pl ]; then
######## echo "translation script is not available, exiting..."
######## exit 1
# if not, get it from GitHub
  git clone https://github.com/geodynamics/specfem2d.wiki.git
  mv -f specfem2d.wiki/translate_user_manual_to_markdown.pl .
  mv -f specfem2d.wiki/*.md .
fi

# runs translation script
./translate_user_manual_to_markdown.pl

echo
echo

