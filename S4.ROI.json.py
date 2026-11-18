import openslide
import cv2
import os
import numpy as np
import json
from tkinter import Tk, filedialog
from pathlib import Path

# ============ 全局变量 ============
drawing = False
start_point = (0, 0)
end_point = (0, 0)
all_rois = []
current_rect = None
current_mouse_pos = (0, 0)  # 记录鼠标实时位置用于放大镜
magnifier_active = True  # 默认开启放大镜


def mouse_callback(event, x, y, flags, param):
    global drawing, start_point, end_point, current_rect, current_mouse_pos

    # 始终更新鼠标位置
    current_mouse_pos = (x, y)

    if event == cv2.EVENT_LBUTTONDOWN:
        drawing = True
        start_point = (x, y)
        end_point = (x, y)
    elif event == cv2.EVENT_MOUSEMOVE:
        if drawing:
            end_point = (x, y)
    elif event == cv2.EVENT_LBUTTONUP:
        drawing = False
        end_point = (x, y)
        x1, y1 = min(start_point[0], end_point[0]), min(start_point[1], end_point[1])
        x2, y2 = max(start_point[0], end_point[0]), max(start_point[1], end_point[1])
        # 过滤掉太小的误触点击
        if (x2 - x1) > 5 and (y2 - y1) > 5:
            current_rect = (x1, y1, x2, y2)
        else:
            current_rect = None


def get_save_directory():
    root = Tk()
    root.withdraw()
    save_dir = filedialog.askdirectory(title='选择保存路径')
    root.destroy()
    return save_dir if save_dir else os.path.join(os.getcwd(), 'roi_results')


def get_screen_resolution():
    """获取屏幕分辨率以适配最大缩略图"""
    root = Tk()
    root.withdraw()
    w, h = root.winfo_screenwidth(), root.winfo_screenheight()
    root.destroy()
    return w, h


def save_roi_data(svs_path, rois, scale_x, scale_y, output_dir, thumb_img):
    if not rois: return

    base_name = Path(svs_path).stem
    save_dir = Path(output_dir) / base_name
    save_dir.mkdir(parents=True, exist_ok=True)

    all_coords = []
    # 复制一份缩略图用于画最终结果
    result_preview = thumb_img.copy()

    for idx, roi_rect in enumerate(rois):
        x1, y1, x2, y2 = roi_rect

        # 映射坐标
        real_x = int(x1 * scale_x)
        real_y = int(y1 * scale_y)
        real_w = int((x2 - x1) * scale_x)
        real_h = int((y2 - y1) * scale_y)

        roi_data = {
            'roi_id': idx + 1,
            'x_start': real_x,
            'y_start': real_y,
            'width': real_w,
            'height': real_h
        }
        all_coords.append(roi_data)

        # 在保存的预览图上绘制带编号的框
        cv2.rectangle(result_preview, (x1, y1), (x2, y2), (0, 255, 0), 2)
        cv2.putText(result_preview, str(idx + 1), (x1, y1 - 5), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0, 255, 0), 2)

    json_path = save_dir / f'{base_name}_coords.json'
    with open(json_path, 'w') as f:
        json.dump(all_coords, f, indent=4)

    cv2.imwrite(str(save_dir / "preview_annotated.png"), result_preview)
    print(f"✅ [保存] {base_name} - {len(rois)} ROIs 已保存")


def annotate_svs(svs_path):
    global start_point, end_point, current_rect, all_rois, drawing, current_mouse_pos, magnifier_active

    # 1. 获取屏幕大小，最大化利用屏幕空间
    screen_w, screen_h = get_screen_resolution()
    max_w = int(screen_w * 0.85)
    max_h = int(screen_h * 0.85)

    try:
        slide = openslide.OpenSlide(svs_path)
    except Exception as e:
        print(f"无法打开文件 {svs_path}: {e}")
        return 'error', [], 1, 1

    # 计算最佳缩略图尺寸
    w_orig, h_orig = slide.dimensions
    aspect_ratio = w_orig / h_orig

    if max_w / aspect_ratio <= max_h:
        thumb_size = (max_w, int(max_w / aspect_ratio))
    else:
        thumb_size = (int(max_h * aspect_ratio), max_h)

    print(f"正在加载 {Path(svs_path).name}，生成高分辨率缩略图 {thumb_size}...")

    thumbnail = slide.get_thumbnail(thumb_size)
    thumbnail_np = np.array(thumbnail.convert('RGB'))
    # OpenCV 使用 BGR，转换颜色空间
    display_image = cv2.cvtColor(thumbnail_np, cv2.COLOR_RGB2BGR)

    # 简单的锐化处理，让组织边缘更清晰
    kernel = np.array([[0, -1, 0], [-1, 5, -1], [0, -1, 0]])
    display_image = cv2.filter2D(display_image, -1, kernel)

    scale_x = w_orig / thumb_size[0]
    scale_y = h_orig / thumb_size[1]

    window_name = f'Main View: {Path(svs_path).name}'
    mag_window_name = 'Magnifier (1:1 Detail)'

    cv2.namedWindow(window_name, cv2.WINDOW_NORMAL)
    cv2.resizeWindow(window_name, thumb_size[0], thumb_size[1])
    cv2.moveWindow(window_name, 50, 50)
    cv2.setMouseCallback(window_name, mouse_callback)

    all_rois = []
    current_rect = None
    drawing = False

    # 提示信息
    print("操作指南:")
    print("  [鼠标左键拖拽] : 框选区域")
    print("  [Space/Enter]  : 确认当前框选")
    print("  [C]            : 清除上一个框选")
    print("  [M]            : 开/关 放大镜")
    print("  [N]            : 下一张图片")
    print("  [Q]            : 退出程序")

    while True:
        # 深拷贝用于绘制UI，不污染原图
        ui_img = display_image.copy()

        # 1. 绘制已保存的ROI (绿色)
        for i, r in enumerate(all_rois):
            cv2.rectangle(ui_img, (r[0], r[1]), (r[2], r[3]), (0, 255, 0), 2)
            # 在框中心画个十字或序号，方便定位
            cx, cy = (r[0] + r[2]) // 2, (r[1] + r[3]) // 2
            cv2.putText(ui_img, str(i + 1), (cx, cy), cv2.FONT_HERSHEY_SIMPLEX, 0.6, (0, 0, 0), 2)

        # 2. 绘制正在拖拽的框 (红色)
        if drawing:
            cv2.rectangle(ui_img, start_point, end_point, (0, 0, 255), 1)
            # 显示实时尺寸
            w_box = abs(end_point[0] - start_point[0])
            h_box = abs(end_point[1] - start_point[1])
            info_txt = f"{w_box}x{h_box}"
            cv2.putText(ui_img, info_txt, (end_point[0] + 5, end_point[1]), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0, 0, 255),
                        1)

        # 3. 绘制十字准星 (辅助对齐)
        mx, my = current_mouse_pos
        cv2.line(ui_img, (mx, 0), (mx, thumb_size[1]), (100, 100, 100), 1)
        cv2.line(ui_img, (0, my), (thumb_size[0], my), (100, 100, 100), 1)

        cv2.imshow(window_name, ui_img)

        # ============ 放大镜逻辑 (核心优化) ============
        if magnifier_active:
            # 映射鼠标坐标到原始大图坐标
            real_mx = int(mx * scale_x)
            real_my = int(my * scale_y)

            # 定义放大镜视窗大小 (例如 512x512 像素)
            mag_w, mag_h = 512, 512

            # 尝试读取区域 (注意边界检查交由openslide处理，但为了速度最好不越界)
            try:
                # 坐标偏移，使鼠标位于放大镜中心
                start_x = max(0, real_mx - mag_w // 2)
                start_y = max(0, real_my - mag_h // 2)

                # 读取 Level 0 (最高分辨率)
                region = slide.read_region((start_x, start_y), 0, (mag_w, mag_h))
                region_np = np.array(region.convert('RGB'))
                region_bgr = cv2.cvtColor(region_np, cv2.COLOR_RGB2BGR)

                # 在放大镜中心画个十字
                cv2.circle(region_bgr, (mag_w // 2, mag_h // 2), 5, (0, 0, 255), 1)
                cv2.line(region_bgr, (mag_w // 2 - 10, mag_h // 2), (mag_w // 2 + 10, mag_h // 2), (0, 0, 255), 1)
                cv2.line(region_bgr, (mag_w // 2, mag_h // 2 - 10), (mag_w // 2, mag_h // 2 + 10), (0, 0, 255), 1)

                cv2.imshow(mag_window_name, region_bgr)
            except Exception:
                pass  # 边缘偶尔读取失败忽略
        else:
            # 如果关闭放大镜，销毁窗口节省资源
            if cv2.getWindowProperty(mag_window_name, cv2.WND_PROP_VISIBLE) >= 1:
                cv2.destroyWindow(mag_window_name)

        # ============ 按键处理 ============
        key = cv2.waitKey(15) & 0xFF
        if key == ord('q'):
            cv2.destroyAllWindows()
            slide.close()
            return 'quit', all_rois, scale_x, scale_y, display_image
        elif key == 13 or key == 32:  # Enter 或 Space 确认选区
            if current_rect:
                all_rois.append(current_rect)
                print(f"已添加 ROI #{len(all_rois)}")
                current_rect = None
        elif key == ord('c'):  # Clear 撤销上一个
            if all_rois:
                all_rois.pop()
                print("已撤销上一个 ROI")
        elif key == ord('m'):  # 切换放大镜
            magnifier_active = not magnifier_active
            print(f"放大镜已{'开启' if magnifier_active else '关闭'}")
        elif key == ord('n'):  # Next
            cv2.destroyAllWindows()
            slide.close()
            return 'next', all_rois, scale_x, scale_y, display_image


def main():
    save_path = get_save_directory()
    print(f"结果将保存至: {save_path}")

    while True:
        root = Tk()
        root.withdraw()
        svs_path = filedialog.askopenfilename(
            title="选择SVS文件",
            filetypes=[('SVS Slide', '*.svs'), ('All Files', '*.*')]
        )
        root.destroy()

        if not svs_path:
            break

        res, rois, sx, sy, thumb = annotate_svs(svs_path)

        # 保存数据
        if rois:
            save_roi_data(svs_path, rois, sx, sy, save_path, thumb)

        if res == 'quit':
            break
        elif res == 'error':
            continue


if __name__ == "__main__":
    main()
