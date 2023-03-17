use std::{error::Error, ptr, slice};

use gdal::{Dataset, DatasetOptions, GdalOpenFlags};
use gdal_sys::GDALRWFlag::GF_Write;

//================================
// RASTER SKELETONIZATION
//================================
// Binary image thinning (skeletonization) in-place.
// Implements Zhang-Suen algorithm.
// http://agcggs680.pbworks.com/f/Zhan-Suen_algorithm.pdf
fn thinning_zs_iteration(
    im: &mut [u8],
    win_x: usize,
    win_y: usize,
    win_w: usize,
    win_h: usize,
    w: usize,
    h: usize,
    iter: i32,
) -> bool {
    let mut diff: bool = false;
    let min_x = if win_x == 0 { 1 } else { win_x };
    let max_x = if win_x + win_w == w {
        w - 1
    } else {
        win_x + win_w
    };
    let min_y = if win_y == 0 { 1 } else { win_y };
    let max_y = if win_y + win_h == h {
        h - 1
    } else {
        win_y + win_h
    };
    for i in min_y..max_y {
        for j in min_x..max_x {
            let p1: u8 = im[i * w + j] & 1;
            let p2: u8 = im[(i - 1) * w + j] & 1;
            let p3: u8 = im[(i - 1) * w + j + 1] & 1;
            let p4: u8 = im[(i) * w + j + 1] & 1;
            let p5: u8 = im[(i + 1) * w + j + 1] & 1;
            let p6: u8 = im[(i + 1) * w + j] & 1;
            let p7: u8 = im[(i + 1) * w + j - 1] & 1;
            let p8: u8 = im[(i) * w + j - 1] & 1;
            let p9: u8 = im[(i - 1) * w + j - 1] & 1;
            let a: u8 = (p2 == 0 && p3 == 1) as u8
                + (p3 == 0 && p4 == 1) as u8
                + (p4 == 0 && p5 == 1) as u8
                + (p5 == 0 && p6 == 1) as u8
                + (p6 == 0 && p7 == 1) as u8
                + (p7 == 0 && p8 == 1) as u8
                + (p8 == 0 && p9 == 1) as u8
                + (p9 == 0 && p2 == 1) as u8;
            let b: u8 = p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9;
            let m1: u8 = if iter == 0 {
                p2 * p4 * p6
            } else {
                p2 * p4 * p8
            };
            let m2: u8 = if iter == 0 {
                p4 * p6 * p8
            } else {
                p2 * p6 * p8
            };
            if a == 1 && (b >= 2 && b <= 6) && m1 == 0 && m2 == 0 && (p1 & 1) == 1 {
                im[i * w + j] |= 2;
                diff = true;
            }
        }
    }

    return diff;
}

fn thinning_zs_post(
    im: &mut [u8],
    win_x: usize,
    win_y: usize,
    win_w: usize,
    win_h: usize,
    w: usize,
) {
    for i in win_y..win_y + win_h {
        for j in win_x..win_x + win_w {
            let marker = im[i * w + j] >> 1;
            let old = im[i * w + j] & 1;
            let new = old & (!marker);
            im[i * w + j] = new;
        }
    }
}

pub fn thinning_zs(im: &mut [u8], w: usize, h: usize) {
    let mut iter = 0;
    let mut diff;
    loop {
        dbg!(iter);
        if dbg!(thinning_zs_iteration(im, 0, 0, w, h, w, h, 0)) {
            thinning_zs_post(im, 0, 0, w, h, w);
            diff = dbg!(thinning_zs_iteration(im, 0, 0, w, h, w, h, 1));
        } else {
            diff = false;
        }
        thinning_zs_post(im, 0, 0, w, h, w);
        if !diff {
            break;
        }
        iter += 1;
    }
}

pub fn thinning_zs_tiled(im: &mut [u8], w: usize, h: usize) {
    let tile_size_x = 256;
    let tile_size_y = 256;
    let ntx = (w + tile_size_x - 1) / tile_size_x;
    let nty = (h + tile_size_y - 1) / tile_size_y;
    let mut tile_done = vec![false; ntx * nty];
    let mut iter = 0;
    loop {
        dbg!(iter);
        let mut diff: bool = false;

        for ti_y in 0..nty {
            for ti_x in 0..ntx {
                if tile_done[ti_y * ntx + ti_x]
                    && (ti_x == 0 || tile_done[ti_y * ntx + ti_x - 1])
                    && (ti_y == 0 || tile_done[(ti_y - 1) * ntx + ti_x])
                    && (ti_x == ntx - 1 || tile_done[ti_y * ntx + ti_x + 1])
                    && (ti_y == nty - 1 || tile_done[(ti_y + 1) * ntx + ti_x])
                {
                    continue;
                }
                // println!("{iter}: {:?}", (ti_x, ti_y));
                let win_x = ti_x * tile_size_x;
                let win_y = ti_y * tile_size_y;
                let win_w = tile_size_x.min(w - win_x);
                let win_h = tile_size_y.min(h - win_y);
                if thinning_zs_iteration(im, win_x, win_y, win_w, win_h, w, h, 0) {
                    diff = true;
                } else {
                    tile_done[ti_y * ntx + ti_x] = true;
                }
            }
        }

        if !diff {
            break;
        }
        for ti_y in 0..nty {
            for ti_x in 0..ntx {
                if tile_done[ti_y * ntx + ti_x] {
                    continue;
                }
                let win_x = ti_x * tile_size_x;
                let win_y = ti_y * tile_size_y;
                let win_w = tile_size_x.min(w - win_x);
                let win_h = tile_size_y.min(h - win_y);
                thinning_zs_post(im, win_x, win_y, win_w, win_h, w);
            }
        }

        // thinning_zs_post(im, 0, 0, w, h, w);
        diff = false;
        for ti_y in 0..nty {
            for ti_x in 0..ntx {
                if tile_done[ti_y * ntx + ti_x]
                    && (ti_x == 0 || tile_done[ti_y * ntx + ti_x - 1])
                    && (ti_y == 0 || tile_done[(ti_y - 1) * ntx + ti_x])
                    && (ti_x == ntx - 1 || tile_done[ti_y * ntx + ti_x + 1])
                    && (ti_y == nty - 1 || tile_done[(ti_y + 1) * ntx + ti_x])
                {
                    continue;
                }
                let win_x = ti_x * tile_size_x;
                let win_y = ti_y * tile_size_y;
                let win_w = tile_size_x.min(w - win_x);
                let win_h = tile_size_y.min(h - win_y);
                if thinning_zs_iteration(im, win_x, win_y, win_w, win_h, w, h, 1) {
                    diff = true;
                } else {
                    tile_done[ti_y * ntx + ti_x] = true;
                }
            }
        }

        if !diff {
            break;
        }
        // thinning_zs_post(im, 0, 0, w, h, w);
        for ti_y in 0..nty {
            for ti_x in 0..ntx {
                if tile_done[ti_y * ntx + ti_x] {
                    continue;
                }
                let win_x = ti_x * tile_size_x;
                let win_y = ti_y * tile_size_y;
                let win_w = tile_size_x.min(w - win_x);
                let win_h = tile_size_y.min(h - win_y);
                thinning_zs_post(im, win_x, win_y, win_w, win_h, w);
            }
        }
        iter += 1;
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    let mut ds = Dataset::open_ex(
        "v.tif",
        DatasetOptions {
            open_flags: GdalOpenFlags::GDAL_OF_UPDATE,
            ..DatasetOptions::default()
        },
    )?;
    let band = ds.rasterband(1)?;
    let mut pixel_space = 0;
    let mut line_space = 0i64;
    let mem = unsafe {
        gdal_sys::GDALGetVirtualMemAuto(
            band.c_rasterband(),
            GF_Write,
            &mut pixel_space as *mut _,
            &mut line_space as *mut _,
            ptr::null::<i8>() as _,
        )
    };
    dbg!(pixel_space);
    dbg!(line_space);
    let (width, height) = band.size();
    assert_eq!(pixel_space, 1);
    assert_eq!(line_space, width as i64);
    let (tile_width, tile_height) = band.block_size();
    dbg!((width, height));
    dbg!((tile_width, tile_height));
    let data = unsafe { gdal_sys::CPLVirtualMemGetAddr(mem) } as *mut u8;
    let len = unsafe { gdal_sys::CPLVirtualMemGetSize(mem) };
    let im = unsafe { slice::from_raw_parts_mut(data, len) };
    // for i in 0..height * width {
    //     if im[i as usize] > 128 {
    //         im[i as usize] = 1
    //     } else {
    //         im[i as usize] = 0
    //     }
    // }

    thinning_zs_tiled(im, width, height);

    unsafe { gdal_sys::CPLVirtualMemFree(mem) };

    ds.flush_cache();

    Ok(())
}
