ffmpeg -framerate 4/1 -i Proj_Combine_%06d.png -c:v libx264 -preset slow -tune animation \
       -pix_fmt yuv420p -s 2400x1200 Proj_Particle_Position.mp4