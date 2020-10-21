import dash_html_components as html

INTRODUCTION_COMPONENT = html.Div(className='pos-f-t col-sm-12', children=[
            html.Div(className='collapse navbar-collapse', id='introduction', children=[
                html.Div(className='col-sm-12 d-flex justify-content-center bg-light p-4', children=[
                    html.Div(className='card', children=[
                        html.Div('Introduction', className='card-header'),
                        html.Div(className='card-body', children=[
                            html.Blockquote(className='blockquote mb-0', children=[
                                html.P('ATAT is a web-based application to study the dynamics of virus '
                                       'reservoir-to-host (R2H) mutation transmissions. The basic requirement of the '
                                       'application is the viral protein sequences of the reservoir and host viruses '
                                       'as an input.'),
                                html.P('The comparative analyses between the co-aligned sequences of the reservoir and '
                                       'host populations is based on a sliding window approach of size nine for '
                                       'statistical significance and data application to the major histocompatibility '
                                       'complex (MHC) and T-cell receptor (TCR) immune response mechanisms.'),
                                html.P('The sequences at each of the aligned overlapping nonamer positions of the '
                                       'respective hosts are classified as four patterns of characteristic diversity '
                                       'motifs, as a basis for quantitative analysis: (i) "index", the most prevalent '
                                       'sequence; (ii) "major" variant, the second most common sequence and the single '
                                       'most prevalent variant of the index, with at least one amino acid '
                                       'mutation; (iii) "minor" variants, multiple different sequences, each with an '
                                       'incidence (percent occurrence) less than that of the major variant; and (iv) '
                                       '"unique" variants, each observed only once.'),
                                html.P('The diversity motifs and their incidences at each of the nonamer positions '
                                       'allow evaluation of the mutation transmission dynamics and selectivity of the '
                                       'sequences in relation to the reservoir or the hosts. This application can be '
                                       'used for a detailed proteome-wide characterization of the composition and '
                                       'incidence of mutations present in the reservoir and host populations for a '
                                       'better understanding of host tropism.'),
                                html.Footer(className='blockquote-footer', children=[
                                    'Dynamics of Influenza A (H5N1) virus protein sequence diversity',
                                    html.Cite(' - Asif M. Khan et al', title='citation')
                                ])

                            ])
                        ])
                    ])
                ])
            ])
        ])