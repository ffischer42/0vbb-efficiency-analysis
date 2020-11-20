#!/bin/bash
screen -dmS calib-worker-1 venv gedet julia calibPrep.jl ../datasets/run0083-run0092-cal-analysis.txt 10 1
screen -dmS calib-worker-2 venv gedet julia calibPrep.jl ../datasets/run0083-run0092-cal-analysis.txt 10 2
screen -dmS calib-worker-3 venv gedet julia calibPrep.jl ../datasets/run0083-run0092-cal-analysis.txt 10 3
screen -dmS calib-worker-4 venv gedet julia calibPrep.jl ../datasets/run0083-run0092-cal-analysis.txt 10 4
screen -dmS calib-worker-5 venv gedet julia calibPrep.jl ../datasets/run0083-run0092-cal-analysis.txt 10 5
screen -dmS calib-worker-6 venv gedet julia calibPrep.jl ../datasets/run0083-run0092-cal-analysis.txt 10 6
screen -dmS calib-worker-7 venv gedet julia calibPrep.jl ../datasets/run0083-run0092-cal-analysis.txt 10 7
screen -dmS calib-worker-8 venv gedet julia calibPrep.jl ../datasets/run0083-run0092-cal-analysis.txt 10 8
screen -dmS calib-worker-9 venv gedet julia calibPrep.jl ../datasets/run0083-run0092-cal-analysis.txt 10 9
screen -dmS calib-worker-10 venv gedet julia calibPrep.jl ../datasets/run0083-run0092-cal-analysis.txt 10 10