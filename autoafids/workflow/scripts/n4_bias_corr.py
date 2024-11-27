# scripts/n4_bias_correction_snakemake.py
import ants


def n4_bias_correction(input_image_path, output_image_path):
    """
    Performs n4 bias field correction.

    Parameters
    ----------
        input_image_path:: str
            Path to the input image

        output_image_path:: str
            Path to the output image

    Returns
    -------
        None
    """
    input_image = ants.image_read(input_image_path)

    # Perform the bias field correction
    output_image = ants.n4_bias_field_correction(input_image)

    ants.image_write(output_image, output_image_path)




def main(snakemake):
    input_image_path = snakemake.input.im
    output_image_path = snakemake.output.corrected_im

    n4_bias_correction(input_image_path, output_image_path)


if __name__ == "__main__":
    main(snakemake)
