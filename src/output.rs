use std::io::ErrorKind;

use image::ImageFormat;

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

pub fn infer_file_format(args: &mut crate::Cli) -> () {
    // LEARN: das & zeichen hier unten war geraten
    //
    if let None = args.image_format {
        if let Some(p) = &args.out.path {
            args.image_format = Some(ImageFormat::from_path(p).unwrap());
        }
    }
}
