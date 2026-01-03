library("ggplot2")
# library("extrafont")

# font_import()
# loadfonts(device="win")

# Journal style definitions (Gyn Onc)

# All units are mm 
oneCol_width = 88
twoCol_width = 190
height_max = 180

# Units are points
font.subscript.size = 5
font.size = 6
font.large.size = 7
font_small_size = font.subscript.size
font_medium_size = font.size
font_large_size = font.large.size

font_family="sans"

font_small = element_text(family = font_family, size=font.subscript.size, color = "black")
font_medium = element_text(family = font_family, size=font.size, color = "black")
font_large =  element_text(family = font_family, size=font.large.size, color = "black")

font_small_md = element_markdown(family = font_family, size=font.subscript.size, color = "black")
font_medium_md = element_markdown(family = font_family, size=font.size, color = "black")
font_large_md =  element_markdown(family = font_family, size=font.large.size, color = "black")

out.format = "pdf"