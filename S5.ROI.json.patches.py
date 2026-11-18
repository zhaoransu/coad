import os
import json
import openslide
from pathlib import Path
from PIL import Image
import time


def crop_patches_from_rois():
    # ================= 配置区域 =================
    # 1. SVS 原始切片文件夹
    svs_dir = Path(r"G:\1.Sprouts\CRC.Pathology\20.WSL\Dataset\High")

    # 2. JSON 坐标文件根目录 (脚本会自动找 子文件夹/xxx_coords.json)
    json_root_dir = Path(r"G:\1.Sprouts\CRC.Pathology\20.WSL\Dataset\High\COAD.High.ROI")

    # 3. 输出目录 (注意：如果G盘满了，请将这里改到其他可用空间大的盘符！)
    output_root = Path(r"G:\1.Sprouts\CRC.Pathology\20.WSL\Dataset\High\COAD.High.ROI.patches")

    # 4. Patch 大小设置
    PATCH_SIZE = 512

    # 是否丢弃边缘不足 224 的残缺块？ (True: 丢弃, False: 也就是只要有图就切，可能小于224)
    # 深度学习通常建议 True，保证输入尺寸一致
    DROP_PARTIAL = True
    # ===========================================

    output_root.mkdir(parents=True, exist_ok=True)

    svs_files = list(svs_dir.glob("*.svs"))
    if not svs_files:
        print(f"❌ 未找到 SVS 文件: {svs_dir}")
        return

    print(f"📂 扫描到 {len(svs_files)} 个样本，准备切割Patches...")
    print("-" * 60)

    total_new_patches_global = 0
    total_skipped_global = 0
    start_time = time.time()

    for svs_path in svs_files:
        slide_name = svs_path.stem
        target_json_path = json_root_dir / slide_name / f"{slide_name}_coords.json"

        if not target_json_path.exists():
            print(f"⚠️  [跳过] {slide_name}: 未找到标注 JSON")
            continue

        print(f"🚀 [处理] {slide_name} ...")

        # 创建该样本的输出目录
        slide_output_dir = output_root / slide_name
        slide_output_dir.mkdir(parents=True, exist_ok=True)

        try:
            with open(target_json_path, 'r', encoding='utf-8') as f:
                rois = json.load(f)

            slide = openslide.OpenSlide(str(svs_path))

            new_patches_in_slide = 0
            skipped_patches_in_slide = 0

            for roi in rois:
                roi_id = roi.get('roi_id', 'unknown')

                # ROI 的全局坐标和尺寸
                roi_x = roi['x_start']
                roi_y = roi['y_start']
                roi_w = roi['width']
                roi_h = roi['height']

                if roi_w < PATCH_SIZE or roi_h < PATCH_SIZE:
                    print(f"    ⚠️ ROI-{roi_id} 尺寸 ({roi_w}x{roi_h}) 小于 {PATCH_SIZE}，无法切割，已跳过。")
                    continue

                # ==============================================================
                # 优化点：在读取整个 ROI 到内存前，先检查一下这个 ROI 的图是不是都已经切完了
                # 如果全切完了，就没必要再把巨大的 ROI 读入内存了，防止无谓的内存消耗
                # ==============================================================
                n_rows = roi_h // PATCH_SIZE
                n_cols = roi_w // PATCH_SIZE

                # 检查此 ROI 下的所有 patch 文件是否都已存在
                all_exist = True
                expected_files = []
                for r in range(n_rows):
                    for c in range(n_cols):
                        global_patch_x = roi_x + c * PATCH_SIZE
                        global_patch_y = roi_y + r * PATCH_SIZE
                        patch_name = f"{slide_name}_ROI{roi_id}_r{r}_c{c}_x{global_patch_x}_y{global_patch_y}.png"
                        save_path = slide_output_dir / patch_name
                        expected_files.append(save_path)
                        if not save_path.exists():
                            all_exist = False

                if all_exist:
                    skipped_patches_in_slide += len(expected_files)
                    continue  # 这个 ROI 已经全切好了，跳过内存读取，直接进行下一个 ROI

                # 如果有没切完的，则将整个 ROI 读入内存
                try:
                    roi_img = slide.read_region((roi_x, roi_y), 0, (roi_w, roi_h))
                    roi_img = roi_img.convert("RGB")
                except Exception as e:
                    print(f"    ❌ 读取 ROI-{roi_id} 失败: {e}")
                    continue

                # === 切割循环 ===
                for r in range(n_rows):
                    for c in range(n_cols):
                        local_x = c * PATCH_SIZE
                        local_y = r * PATCH_SIZE

                        global_patch_x = roi_x + local_x
                        global_patch_y = roi_y + local_y

                        patch_name = (
                            f"{slide_name}_ROI{roi_id}_"
                            f"r{r}_c{c}_"
                            f"x{global_patch_x}_y{global_patch_y}.png"
                        )
                        save_path = slide_output_dir / patch_name

                        # 🌟 核心断点续传逻辑：如果已经存在，就跳过不保存
                        if save_path.exists():
                            skipped_patches_in_slide += 1
                            continue

                        # 裁剪并保存
                        box = (local_x, local_y, local_x + PATCH_SIZE, local_y + PATCH_SIZE)
                        patch = roi_img.crop(box)

                        try:
                            patch.save(save_path)
                            new_patches_in_slide += 1
                        except OSError as e:
                            # 捕获磁盘满的错误，抛出到外层停止运行
                            raise Exception(f"保存图片失败，可能是硬盘满了: {e}")

            slide.close()
            total_new_patches_global += new_patches_in_slide
            total_skipped_global += skipped_patches_in_slide

            if skipped_patches_in_slide > 0 and new_patches_in_slide == 0:
                print(f"    ⏭️  已跳过 {skipped_patches_in_slide} 个已存在的 Patches")
            else:
                print(f"    ✅ 新生成 {new_patches_in_slide} 个 | 跳过 {skipped_patches_in_slide} 个")

        except Exception as e:
            print(f"    ❌ 处理出错: {e}")
            if "硬盘满了" in str(e) or "No space left" in str(e):
                print("\n🚨 严重错误：磁盘空间不足！程序安全中止以防止数据损坏。请清理硬盘空间后再次运行本脚本，它会自动从当前位置继续。")
                break  # 磁盘满了就直接中断大循环

    duration = time.time() - start_time
    print("-" * 60)
    print(f"🎉 任务结束！")
    print(f"📊 本次新产出: {total_new_patches_global}")
    print(f"⏭️ 本次共跳过: {total_skipped_global}")
    print(f"⏱️ 总耗时: {duration:.2f} 秒")
    print(f"📁 输出位置: {output_root}")


if __name__ == "__main__":
    # os.add_dll_directory(r"C:\openslide-win64\bin") # 如果报错 DLL load failed 请取消注释
    crop_patches_from_rois()