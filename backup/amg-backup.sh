#!/bin/bash
tar -cvf myamg-backup-$(date +%Y-%m-%d-%H-%M-%S).tar.gz \
			  ../include/ \
                          ../src/*.c ../src/Makefile \
			  ../util/ \
			  ../backup/*.sh \
			  ../make.inc \
			  ../test/ \
			  ../dat/ \
			  ../output/ \
			  ../example/ \
			  ../lib/ \
			  ../doc/ \
			  --exclude=*~ \
			  --exclude=*.swp \
			  --exclude=../lib/* \
			  --exclude=../output/* \
			  --exclude=../dat/*
