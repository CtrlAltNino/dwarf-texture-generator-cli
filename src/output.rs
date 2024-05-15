use std::{
    io::{self, Cursor, ErrorKind, Write},
    path::Path,
};

use image::{ImageError, ImageFormat, ImageResult};

pub fn parse_file_format(format_ext: &str) -> Result<ImageFormat, std::io::Error> {
    //LEARN: is das richtig so, let fail?, am ende to_string?, Wie wuerde man es machen, wenn man
    //fail nach dem let aber for dem Err nochmal aendert?
    let fail;
    if let Some(iform) = ImageFormat::from_extension(format_ext) {
        if iform.writing_enabled() {
            return Ok(iform);
        } else {
            fail = "Cannot encode this image format. Try:";
        }
    } else {
        fail = "Image format not recognized. Try:";
    }
    Err(std::io::Error::new(
        ErrorKind::Other,
        ImageFormat::all()
            .filter(ImageFormat::can_write)
            //LEARN: to_owned vs to_string here
            .fold(fail.to_string(), |accum, iform| {
                accum + " " + iform.extensions_str().first().unwrap()
            }),
    ))
}

pub fn infer_file_format(args: &mut crate::Cli) -> Result<(), ImageError> {
    // LEARN: das & zeichen hier unten war geraten
    if let None = args.image_format {
        if let Some(p) = &args.out.path {
            args.image_format = Some(ImageFormat::from_path(p)?);
        } else {
            args.image_format = Some(ImageFormat::Png);
        }
    }
    Ok(())
}

//TODO: implement error handling here?
pub fn save_to_file(
    image: image::ImageBuffer<image::Rgb<u8>, Vec<u8>>,
    iformat: ImageFormat,
    pth: &Path,
) -> () {
    image.save_with_format(pth, iformat);
}

pub fn spit_to_stdout(
    image: image::ImageBuffer<image::Rgb<u8>, Vec<u8>>,
    iformat: ImageFormat,
) -> () {
    let mut cursor = Cursor::new(Vec::new());

    image.write_to(&mut cursor, iformat).unwrap();
    //TODO gracefully terminate here
    io::stdout().write_all(&cursor.into_inner()).unwrap();
    io::stdout().flush().unwrap();
}
