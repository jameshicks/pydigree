#!/bin/bash

nosetests -v \
	  --with-coverage \
	  --cover-inclusive \
	  --cover-html \
	  --cover-package=pydigree \
	  --cover-branches
