#!/bin/sh

raycreate room 1
cp room.ply room2.ply
raydecimate room2.ply 10 cm
raytranslate room2_decimated.ply 1,2,3
rayrotate room2_decimated.ply 0,0,-50
rayrestore room2_decimated.ply 10 cm room.ply

