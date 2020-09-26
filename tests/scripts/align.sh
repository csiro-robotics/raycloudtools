#!/bin/sh

raycreate room 1
cp room.ply room2.ply
raytranslate room2.ply 1,2,3
rayrotate room2.ply 0,0,35
rayalign room.ply room2.ply

