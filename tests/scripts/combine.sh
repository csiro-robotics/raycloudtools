#!/bin/sh

raycreate room 1
cp room.ply room2.ply
raytranslate room2.ply 0,0,1
rayrotate room2.ply 0,0,35
raycombine min room.ply room2.ply 1 rays
