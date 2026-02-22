library(magick)
library(hexSticker)

# 1. 读入并处理简笔画图案
# 关键修复：这里的白底会撑出六边形（因为 hexSticker 的 s_width 如果缩放没控制好边缘的白色实底）
# 我们彻底把它变成一头"没有底色、完全透明"的图，让六边形(h_fill)去处理底色
img_raw <- image_read("/Users/leoarrow/.gemini/antigravity/brain/610eaede-f728-455b-b2b1-2cece78f539e/simple_doodle_lion_1771730123312.png")
img_clean <- image_trim(img_raw)
img_trans <- image_transparent(img_clean, "white", fuzz = 8)

img_ready <- image_trim(img_trans)

image_write(img_ready, "man/figures/icon_processed.png")

# 2. 生成最终 hexSticker
# 确保六边形背景和底层透明能够融合
sticker(
  subplot = "man/figures/icon_processed.png", package = "leo.ukb", 
  p_size = 22, p_y = 1.4, p_color = "#336699", 
  s_x = 1, s_y = 0.75, s_width = 0.50, s_height = 0.50,  # 进一步缩小了芯的大小，绝对不会凸出来
  h_fill = "#FFFFFF", h_color = "#336699", 
  filename = "man/figures/logo.png"
)
