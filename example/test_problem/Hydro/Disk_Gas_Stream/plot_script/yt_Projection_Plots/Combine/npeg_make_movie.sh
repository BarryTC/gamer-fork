ffmpeg -framerate 2/1 -i Combine_%06d.png -c:v libx264 -preset slow -tune film \
       -pix_fmt yuv420p -s 2400x2400 Proj_Star_Gas.mp4