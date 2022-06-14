use bgzip::read::BGzReader;
use std::io::{stdout, BufRead, BufReader, BufWriter, Write};
use std::{env, fs};

fn main() {
    eprintln!("[usage]: vcf2gt <vcf.gz>\n");

    let args: Vec<String> = env::args().collect();
    let vcf_name = check_args(&args);

    let vcf =
        BGzReader::new(fs::File::open(vcf_name).unwrap()).expect("Your vcf couldn't be opened!");

    let out = stdout();
    let mut out = BufWriter::new(out.lock());

    for line in BufReader::new(vcf).lines() {
        let line = line.unwrap();
        if line.starts_with("#") {
            continue;
        } else {
            let cols: Vec<&str> = line.split('\t').collect();
            write!(out, "{}", &cols[0]).unwrap();
            write!(out, " {}", &cols[1]).unwrap();

            let gts = cols[9..].iter()
                               .map(|col| col.split(':').collect::<Vec<&str>>()[0])
                               .map(|gt| translate_gt(gt))
                               .collect::<Vec<i32>>();
            for gt in gts {
                write!(out, " {}", gt).unwrap();
            }
            write!(out, "\n").unwrap();
        }
    }
}

fn check_args(args: &[String]) -> &str {
    if args.len() < 2 {
        panic!("Please input your VCF!");
    } else {
        &args[1]
    }
}

fn translate_gt(gt: &str) -> i32 {
    if gt == "0/0" {
        0
    } else if gt == "0/1" {
        1
    } else if gt == "1/1" {
        2
    } else if gt == "./." {
        9
    } else {
        panic!("{} coudn't be recognized!!", gt);
    }
}
